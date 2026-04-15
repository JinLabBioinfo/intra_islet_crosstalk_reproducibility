# Author: Jiachen Sun

source("/home/jxs2269/islet_cosmx/for_jasmine/src/pipeline_reproducibility_common/_utility_funs.r")

# Pre-process part 1: data loading, quality control, sketching, clustering, and integration projection.
obj_ak090t <- pre_process_islet_part_1("AK090T")
obj_al189t <- pre_process_islet_part_1("AL189T")
obj_am010 <- pre_process_islet_part_1("AM010")

# Pre-process part 2: coarse cell typing and c3 graph construction.
obj_ak090t <- pre_process_islet_part_2(obj_ak090t, "AK090T")
obj_al189t <- pre_process_islet_part_2(obj_al189t, "AL189T")
obj_am010 <- pre_process_islet_part_2(obj_am010, "AM010")

# Pre-process part 3: manual selection of exocrine regions.
pre_process_islet_part_3(obj_ak090t, "AK090T", c(443, 442, 421, 530, 513, 514, 562, 563, 193, 191, 210, 311, 297, 298, 18, 19, 122, 121))
pre_process_islet_part_3(obj_al189t, "AL189T", c(39, 144, 159, 280, 310, 369, 450, 184, 273, 530, 455, 480, 565, 205, 488, 489, 350, 659))
pre_process_islet_part_3(obj_am010, "AM010", c(25, 37, 38, 113, 149, 198, 200, 202, 233, 304, 332, 365, 367, 390, 394, 399, 439, 443))
