//
// Created by 63175 on 2019/9/10.
//

#ifndef BWAMEZ_DEBUG_0_5_4_BWAMEZ_H
#define BWAMEZ_DEBUG_0_5_4_BWAMEZ_H

#include "bwa/bwamem.h"

#ifdef __cplusplus
extern "C" {
#endif

	void mez_process_seqs(const mem_opt_t *opt, const bwt_t *bwt, const bntseq_t *bns, const uint8_t *pac,
			int64_t n_processed, int n, bseq1_t *seqs, const mem_pestat_t *pes0);

#ifdef __cplusplus
}
#endif

#endif //BWAMEZ_DEBUG_0_5_4_BWAMEZ_H
