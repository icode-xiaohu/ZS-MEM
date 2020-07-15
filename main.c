#include <zlib.h>
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <ctype.h>
#include <math.h>
#include <getopt.h>

#include "bwa/kseq.h"
#include "bwa/bwa.h"
#include "bwa/bwamem.h"
#include "bwa/kvec.h"
#include "bwa/utils.h"
#include "bwa/bntseq.h"
#include "bwamez.h"
#include "time_prof.h"

KSEQ_DECLARE(gzFile)

extern unsigned char nst_nt4_table[256];

void *kopen(const char *fn, int *_fd);
int kclose(void *a);
void kt_pipeline(int n_threads, void *(*func)(void*, int, void*), void *shared_data, int n_steps);

double tprof[256][4];

typedef struct {
	kseq_t *ks, *ks2;
	mem_opt_t *opt;
	mem_pestat_t *pes0;
	int64_t n_processed;
	int copy_comment, actual_chunk_size;
	bwaidx_t *idx;
} ktp_aux_t;

typedef struct {
	ktp_aux_t *aux;
	int n_seqs;
	bseq1_t *seqs;
} ktp_data_t;

static void *process(void *shared, int step, void *_data)
{
	ktp_aux_t *aux = (ktp_aux_t*)shared;
	ktp_data_t *data = (ktp_data_t*)_data;
	int i;
	if (step == 0) {
		ktp_data_t *ret;
		int64_t size = 0;
		ret = calloc(1, sizeof(ktp_data_t));
		ret->seqs = bseq_read(aux->actual_chunk_size, &ret->n_seqs, aux->ks, aux->ks2);
		if (ret->seqs == 0) {
			free(ret);
			return 0;
		}
		if (!aux->copy_comment)
			for (i = 0; i < ret->n_seqs; ++i) {
				free(ret->seqs[i].comment);
				ret->seqs[i].comment = 0;
			}
		for (i = 0; i < ret->n_seqs; ++i) size += ret->seqs[i].l_seq;
		if (bwa_verbose >= 3)
			fprintf(stderr, "[M::%s] read %d sequences (%ld bp), processed %ld sequences...\n", __func__, ret->n_seqs, (long)size, aux->n_processed);
		return ret;
	} else if (step == 1) {
		const mem_opt_t *opt = aux->opt;
		const bwaidx_t *idx = aux->idx;
		if (opt->flag & MEM_F_SMARTPE) {
			bseq1_t *sep[2];
			int n_sep[2];
			mem_opt_t tmp_opt = *opt;
			bseq_classify(data->n_seqs, data->seqs, n_sep, sep);
			if (bwa_verbose >= 3)
				fprintf(stderr, "[M::%s] %d single-end sequences; %d paired-end sequences\n", __func__, n_sep[0], n_sep[1]);
			if (n_sep[0]) {
				tmp_opt.flag &= ~MEM_F_PE;
				mem_process_seqs(&tmp_opt, idx->bwt, idx->bns, idx->pac, aux->n_processed, n_sep[0], sep[0], 0);
				for (i = 0; i < n_sep[0]; ++i)
					data->seqs[sep[0][i].id].sam = sep[0][i].sam;
			}
			if (n_sep[1]) {
				tmp_opt.flag |= MEM_F_PE;
				mem_process_seqs(&tmp_opt, idx->bwt, idx->bns, idx->pac, aux->n_processed + n_sep[0], n_sep[1], sep[1], aux->pes0);
				for (i = 0; i < n_sep[1]; ++i)
					data->seqs[sep[1][i].id].sam = sep[1][i].sam;
			}
			free(sep[0]); free(sep[1]);
		} else {
			mez_process_seqs(opt, idx->bwt, idx->bns, idx->pac, aux->n_processed, data->n_seqs, data->seqs, aux->pes0);
		}
		aux->n_processed += data->n_seqs;
		return data;
	} else if (step == 2) {
		for (i = 0; i < data->n_seqs; ++i) {
			if (data->seqs[i].sam) err_fputs(data->seqs[i].sam, stdout);
			free(data->seqs[i].name); free(data->seqs[i].comment);
			free(data->seqs[i].seq); free(data->seqs[i].qual);
			free(data->seqs[i].sam);
		}
		free(data->seqs); free(data);
		return 0;
	}
	return 0;
}

static int update_scoring(const char *mode, mem_opt_t *opt, const mem_opt_t *opt0)
{
	if (mode) {
		if (strcmp(mode, "intractg") == 0) {
			if (!opt0->o_del) opt->o_del = 16;
			if (!opt0->o_ins) opt->o_ins = 16;
			if (!opt0->b) opt->b = 9;
			if (!opt0->pen_clip5) opt->pen_clip5 = 5;
			if (!opt0->pen_clip3) opt->pen_clip3 = 5;
		} else if (strcmp(mode, "pacbio") == 0 || strcmp(mode, "pbref") == 0 || strcmp(mode, "ont2d") == 0) {
			if (!opt0->o_del) opt->o_del = 1;
			if (!opt0->e_del) opt->e_del = 1;
			if (!opt0->o_ins) opt->o_ins = 1;
			if (!opt0->e_ins) opt->e_ins = 1;
			if (!opt0->b) opt->b = 1;
			if (opt0->split_factor == 0.) opt->split_factor = 10.;
			if (strcmp(mode, "ont2d") == 0) {
				if (!opt0->min_chain_weight) opt->min_chain_weight = 20;
				if (!opt0->min_seed_len) opt->min_seed_len = 14;
				if (!opt0->pen_clip5) opt->pen_clip5 = 0;
				if (!opt0->pen_clip3) opt->pen_clip3 = 0;
			} else {
				if (!opt0->min_chain_weight) opt->min_chain_weight = 40;
				if (!opt0->min_seed_len) opt->min_seed_len = 17;
				if (!opt0->pen_clip5) opt->pen_clip5 = 0;
				if (!opt0->pen_clip3) opt->pen_clip3 = 0;
			}
		} else {
			fprintf(stderr, "[E::%s] unknown read type '%s'\n", __func__, mode);
			return 1;
		}
	} else {
		if (opt0->a) { // matching score is changed
			if (!opt0->b) opt->b *= opt->a;
			if (!opt0->T) opt->T *= opt->a;
			if (!opt0->o_del) opt->o_del *= opt->a;
			if (!opt0->e_del) opt->e_del *= opt->a;
			if (!opt0->o_ins) opt->o_ins *= opt->a;
			if (!opt0->e_ins) opt->e_ins *= opt->a;
			if (!opt0->zdrop) opt->zdrop *= opt->a;
			if (!opt0->pen_clip5) opt->pen_clip5 *= opt->a;
			if (!opt0->pen_clip3) opt->pen_clip3 *= opt->a;
			if (!opt0->pen_unpaired) opt->pen_unpaired *= opt->a;
		}
	}
	return 0;
}

static mem_opt_t *opt_init() {
	mem_opt_t *opt = mem_opt_init();
	opt->min_reads = 2;
	opt->max_mis = 0.05;
	return opt;
}

static int usage(mem_opt_t *opt) {
	fprintf(stderr, "\n");
	fprintf(stderr, "Usage: zs-mem [options] <idxbase> <in1.fq>\n\n");
	fprintf(stderr, "tune options:\n\n");
	fprintf(stderr, "Zip-seeding options:\n\n");
	fprintf(stderr, "    --min-reads INT    minimum reads to form a Consensus sequence (CS) [%d]\n", opt->min_reads);
	fprintf(stderr, "    --max-mis   FLOAT  maximum mismatches between read and CS [%.2f]\n", opt->max_mis);
	fprintf(stderr, "    -t          INT    number of threads [%d]\n", opt->n_threads);
	fprintf(stderr, "    -k          INT    minimum seed length [%d]\n", opt->min_seed_len);
	fprintf(stderr, "    -r          FLOAT  look for internal seeds inside a seed longer than {-k} * FLOAT [%g]\n", opt->split_factor);
	fprintf(stderr, "    -y          INT    seed occurrence for the 3rd round seeding [%ld]\n", (long)opt->max_mem_intv);
	fprintf(stderr, "    -c          INT    skip seeds with more than INT occurrences [%d]\n", opt->max_occ);

	fprintf(stderr, "\n");
	fprintf(stderr, "Mapping options:\n\n");
	fprintf(stderr, "       -w INT        band width for banded alignment [%d]\n", opt->w);
	fprintf(stderr, "       -d INT        off-diagonal X-dropoff [%d]\n", opt->zdrop);
	fprintf(stderr, "       -D FLOAT      drop chains shorter than FLOAT fraction of the longest overlapping chain [%.2f]\n", opt->drop_ratio);
	fprintf(stderr, "       -W INT        discard a chain if seeded bases shorter than INT [0]\n");
//	fprintf(stderr, "       -m INT        perform at most INT rounds of mate rescues for each read [%d]\n", opt->max_matesw);
//	fprintf(stderr, "       -S            skip mate rescue\n");
//	fprintf(stderr, "       -P            skip pairing; mate rescue performed unless -S also in use\n");

	fprintf(stderr, "\n");
	fprintf(stderr, "Scoring options:\n\n");
	fprintf(stderr, "       -A INT        score for a sequence match, which scales options -TdBOELU unless overridden [%d]\n", opt->a);
	fprintf(stderr, "       -B INT        penalty for a mismatch [%d]\n", opt->b);
	fprintf(stderr, "       -O INT[,INT]  gap open penalties for deletions and insertions [%d,%d]\n", opt->o_del, opt->o_ins);
	fprintf(stderr, "       -E INT[,INT]  gap extension penalty; a gap of size k cost '{-O} + {-E}*k' [%d,%d]\n", opt->e_del, opt->e_ins);
	fprintf(stderr, "       -L INT[,INT]  penalty for 5'- and 3'-end clipping [%d,%d]\n", opt->pen_clip5, opt->pen_clip3);
//	fprintf(stderr, "       -U INT        penalty for an unpaired read pair [%d]\n\n", opt->pen_unpaired);
	fprintf(stderr, "       -x STR        read type. Setting -x changes multiple parameters unless overridden [null]\n");
	fprintf(stderr, "                     pacbio: -k17 -W40 -r10 -A1 -B1 -O1 -E1 -L0  (PacBio reads to ref)\n");
	fprintf(stderr, "                     ont2d: -k14 -W20 -r10 -A1 -B1 -O1 -E1 -L0  (Oxford Nanopore 2D-reads to ref)\n");
	fprintf(stderr, "                     intractg: -B9 -O16 -L5  (intra-species contigs to ref)\n");

	fprintf(stderr, "\n");
	fprintf(stderr, "I/O options:\n\n");
//	fprintf(stderr, "       -p            smart pairing (ignoring in2.fq)\n");
	fprintf(stderr, "       -R STR        read group header line such as '@RG\\tID:foo\\tSM:bar' [null]\n");
	fprintf(stderr, "       -H STR/FILE   insert STR to header if it starts with @; or insert lines in FILE [null]\n");
	fprintf(stderr, "       -o FILE       sam file to output results to [stdout]\n");
	fprintf(stderr, "       -j            treat ALT contigs as part of the primary assembly (i.e. ignore <idxbase>.alt file)\n");
	fprintf(stderr, "       -5            for split alignment, take the alignment with the smallest coordinate as primary\n");
	fprintf(stderr, "       -q            don't modify mapQ of supplementary alignments\n");
	fprintf(stderr, "       -K INT        process INT input bases in each batch regardless of nThreads (for reproducibility) []\n");
	fprintf(stderr, "       -v INT        verbosity level: 1=error, 2=warning, 3=message, 4+=debugging [%d]\n", bwa_verbose);
	fprintf(stderr, "       -T INT        minimum score to output [%d]\n", opt->T);
	fprintf(stderr, "       -h INT[,INT]  if there are <INT hits with score >80%% of the max score, output all in XA [%d,%d]\n",
			opt->max_XA_hits, opt->max_XA_hits_alt);
	fprintf(stderr, "       -a            output all alignments for SE or unpaired PE\n");
	fprintf(stderr, "       -C            append FASTA/FASTQ comment to SAM output\n");
	fprintf(stderr, "       -V            output the reference FASTA header in the XR tag\n");
	fprintf(stderr, "       -Y            use soft clipping for supplementary alignments\n");
	fprintf(stderr, "       -M            mark shorter split hits as secondary\n");
	fprintf(stderr, "\n");
	return 1;
}

int main(int argc, char **argv) {
	double rtime = realtime();
	mem_opt_t *opt, opt0;
	int fd, i, c, lo_index, ignore_alt = 0, no_mt_io = 0;
	int fixed_chunk_size = -1;
	gzFile fp = 0;
	char *p;
	char *rg_line = 0, *hdr_line = 0;
	char *mode = 0;
	void *ko = 0;
	mem_pestat_t pes[4];
	ktp_aux_t aux;

	memset(&aux, 0, sizeof(ktp_aux_t));
	memset(pes, 0, 4 * sizeof(mem_pestat_t));
	for (i = 0; i < 4; ++i) pes[i].failed = 1;

	aux.opt = opt = opt_init();
	memset(&opt0, 0, sizeof(mem_opt_t));
	const char* short_opts = "51qpaMCSPVYjk:c:v:s:r:t:R:A:B:O:E:U:w:L:d:T:Q:D:m:I:N:o:f:W:x:G:h:y:K:X:H:";
	const struct option long_opts[] = {
			{"min-reads", required_argument, NULL, 0},
			{"max-mis", required_argument, NULL, 0},
			{NULL, 0, NULL, 0}
	};
	while((c=getopt_long(argc, argv, short_opts, long_opts, &lo_index)) >= 0) {

		/* Zip-seding options */
		if(c == 0) {
			if(!strcmp(long_opts[lo_index].name, "min-reads")) {
				opt->min_reads = atoi(optarg);
			}
			else if(!strcmp(long_opts[lo_index].name, "max-mis")) {
				opt->max_mis = atof(optarg);
			}
		}
		else if (c == 't') opt->n_threads = atoi(optarg), opt->n_threads = opt->n_threads > 1 ?opt->n_threads :1;
		else if (c == 'k') opt->min_seed_len = atoi(optarg);
		else if (c == 'r') opt->split_factor = atof(optarg);
		else if (c == 'y') opt->max_mem_intv = atol(optarg);
		else if (c == 'c') opt->max_occ = atoi(optarg);

		/* Mapping options */
		else if (c == 'w') opt->w = atoi(optarg);
		else if (c == 'd') opt->zdrop = atoi(optarg);
		else if (c == 'D') opt->drop_ratio = atof(optarg);
		else if (c == 'W') opt->min_chain_weight = atoi(optarg);
		/* Scoring options */
		else if (c == 'A') opt->a = atoi(optarg), opt0.a = 1;
		else if (c == 'B') opt->b = atoi(optarg), opt0.b = 1;
		else if (c == 'O') {
			opt->o_del = opt->o_ins = strtol(optarg, &p, 10);
			if (*p != 0 && ispunct(*p) && isdigit(p[1]))
				opt->o_ins = strtol(p+1, &p, 10);
			opt0.o_del = opt0.o_ins = 1;
		}
		else if (c == 'E') {
			opt->e_del = opt->e_ins = strtol(optarg, &p, 10);
			if (*p != 0 && ispunct(*p) && isdigit(p[1]))
				opt->e_ins = strtol(p+1, &p, 10);
			opt0.e_del = opt0.e_ins = 1;
		}
		else if (c == 'L') {
			opt->pen_clip5 = opt->pen_clip3 = strtol(optarg, &p, 10);
			if (*p != 0 && ispunct(*p) && isdigit(p[1]))
				opt->pen_clip3 = strtol(p+1, &p, 10);
			opt0.pen_clip5 = opt0.pen_clip3 = 1;
		}
		else if (c == 'x') mode = optarg;

		/* I/O options */
		else if (c == 'R') {
			if ((rg_line = bwa_set_rg(optarg)) == 0) return 1; // FIXME: memory leak
		}
		else if (c == 'H') {
			if (optarg[0] != '@') {
				FILE *fp;
				if ((fp = fopen(optarg, "r")) != 0) {
					char *buf;
					buf = calloc(1, 0x10000);
					while (fgets(buf, 0xffff, fp)) {
						i = strlen(buf);
						assert(buf[i-1] == '\n'); // a long line
						buf[i-1] = 0;
						hdr_line = bwa_insert_header(buf, hdr_line);
					}
					free(buf);
					fclose(fp);
				}
			} else hdr_line = bwa_insert_header(optarg, hdr_line);
		}
		else if (c == 'o' || c == 'f') xreopen(optarg, "wb", stdout);
		else if (c == 'j') ignore_alt = 1;
		else if (c == '5') opt->flag |= MEM_F_PRIMARY5 | MEM_F_KEEP_SUPP_MAPQ; // always apply MEM_F_KEEP_SUPP_MAPQ with -5
		else if (c == 'q') opt->flag |= MEM_F_KEEP_SUPP_MAPQ;
		else if (c == 'K') fixed_chunk_size = atoi(optarg);
		else if (c == 'v') bwa_verbose = atoi(optarg);
		else if (c == 'T') opt->T = atoi(optarg), opt0.T = 1;
		else if (c == 'h') {
			opt->max_XA_hits = opt->max_XA_hits_alt = strtol(optarg, &p, 10);
			if (*p != 0 && ispunct(*p) && isdigit(p[1]))
				opt->max_XA_hits_alt = (int)strtol(p+1, &p, 10);
		}
		else if (c == 'a') opt->flag |= MEM_F_ALL;
		else if (c == 'C') aux.copy_comment = 1;
		else if (c == 'V') opt->flag |= MEM_F_REF_HDR;
		else if (c == 'Y') opt->flag |= MEM_F_SOFTCLIP;
		else if (c == 'M') opt->flag |= MEM_F_NO_MULTI;

		/* For debugging */
		else if (c == '1') no_mt_io = 1;

		else {
			usage(opt); free(opt);
			return 1;
		}
	}

	/* Opt[i] and Opt[i+1] should be ref and reads */
	if (optind + 1 >= argc || optind + 2 < argc) {
		usage(opt); free(opt);
		return 1;
	}

	/* Insert header */
	if (rg_line) {
		hdr_line = bwa_insert_header(rg_line, hdr_line);
		free(rg_line);
	}

	/* Update SW scoring */
	int mode_ret = update_scoring(mode, opt, &opt0);
	if(mode_ret) {
		free(opt);
		return 1;
	}
	bwa_fill_scmat(opt->a, opt->b, opt->mat);

	/* Load ref-index */
	aux.idx = bwa_idx_load_from_shm(argv[optind]);
	if (aux.idx == 0) {
		if ((aux.idx = bwa_idx_load(argv[optind], BWA_IDX_ALL)) == 0) return 1; // FIXME: memory leak
	} else if (bwa_verbose >= 3)
		fprintf(stderr, "[M::%s] load the bwa index from shared memory\n", __func__);
	if (ignore_alt) {
		for (i = 0; i < aux.idx->bns->n_seqs; ++i)
			aux.idx->bns->anns[i].is_alt = 0;
	}

	/* Input reads */
	aux.actual_chunk_size = fixed_chunk_size > 0? fixed_chunk_size : opt->chunk_size * opt->n_threads;
	ko = kopen(argv[optind + 1], &fd);
	if (ko == 0) {
		if (bwa_verbose >= 1) fprintf(stderr, "[E::%s] fail to open file `%s'.\n", __func__, argv[optind + 1]);
		return 1;
	}
	fp = gzdopen(fd, "r");
	aux.ks = kseq_init(fp);

	/* Output SAM header */
	bwa_print_sam_hdr(aux.idx->bns, hdr_line);

	/* Start mapping */
	memset(tprof, 0, sizeof(tprof));
	tprof[TP_Z_MAPPING][0] = cputime(); tprof[TP_Z_MAPPING][1] = realtime();
	kt_pipeline(no_mt_io ? 1 : 2, process, &aux, 3);
	tprof[TP_Z_MAPPING][2] = cputime() - tprof[TP_Z_MAPPING][0]; tprof[TP_Z_MAPPING][3] = realtime() - tprof[TP_Z_MAPPING][1];

	/* Free allocated memory */
	free(hdr_line);
	free(opt);
	bwa_idx_destroy(aux.idx);
	kseq_destroy(aux.ks);
	err_gzclose(fp); kclose(ko);

	fprintf(stderr, "Time profiling\n");
	fprintf(stderr, "    Overall mapping: %.2f CPU sec; %.2f real sec\n", tprof[TP_Z_MAPPING][2], tprof[TP_Z_MAPPING][3]);
	fprintf(stderr, "    Seeding:         %.2f CPU sec; %.2f real sec\n", tprof[TP_Z_SEEDING][2], tprof[TP_Z_SEEDING][3]);
	fprintf(stderr, "    Non-seeding:     %.2f CPU sec; %.2f real sec\n", tprof[TP_Z_DP][2], tprof[TP_Z_DP][3]);
	fprintf(stderr, "[ZS-MEM: %s] CPU time: %.2f sec; real time: %.2f sec\n",
			__func__, cputime(), realtime()-rtime);
	return 0;
}