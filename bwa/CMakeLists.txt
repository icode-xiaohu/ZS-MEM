cmake_minimum_required(VERSION 3.13)
project(zip_seeding_debug_5_8 C)

set(CMAKE_C_STANDARD 11)

include_directories(.)

add_executable(zip_seeding_debug_5_8
		bamlite.c
		bamlite.h
		bntseq.c
		bntseq.h
		bwa.c
		bwa.h
		bwamem.c
		bwamem.h
		bwamem_extra.c
		bwamem_pair.c
		bwamez.c
		bwamez.h
		bwape.c
		bwase.c
		bwase.h
		bwaseqio.c
		bwashm.c
		bwt.c
		bwt.h
		bwt_gen.c
		bwt_lite.c
		bwt_lite.h
		bwtaln.c
		bwtaln.h
		bwtgap.c
		bwtgap.h
		bwtindex.c
		bwtsw2.h
		bwtsw2_aux.c
		bwtsw2_chain.c
		bwtsw2_core.c
		bwtsw2_main.c
		bwtsw2_pair.c
		example.c
		fastmap.c
		is.c
		kbtree.h
		khash.h
		kopen.c
		kseq.h
		ksort.h
		kstring.c
		kstring.h
		ksw.c
		ksw.h
		kthread.c
		kvec.h
		main.c
		malloc_wrap.c
		malloc_wrap.h
		maxk.c
		pemerge.c
		QSufSort.c
		QSufSort.h
		rle.c
		rle.h
		rope.c
		rope.h
		utils.c
		utils.h
		zips_main.c)
