# ZENODO ##########

# TODO:

# for preproc: aim for 1 directory per series
#    case of too many tif files: move them to another directory "toCurate" (already sorted by Theo)
#    case of csv only (have been copied from the `preprocess` directory - write a script that find all export files 
# recursively in a directory, and copy the parent subdirectory to the curated subfolder)
# then rerun archiving

# for d in series_dir:
#   mkdir d in output_dir
#   for f in tree(d)
#     ln -s f output_dir/d
#     for f2 in basenmae with tif but no csv
#       rm f2

# for d in series_dir:
#   mkdir d in output_dir
#   for f in d:
#     p = find f in preproc[d]
#     ln -s dirname(p) output_dir/d



(function(export_dir = "/scicore/home/nimwegen/GROUP/MM_Data/Thomas/_GCSlacArticle/Julou_2022_GCS_RawImages",
          .force=FALSE) {
# this anonymous function exports raw data (one dataset per compressed file) after the following scheme:
# - traverse R data to list all datasets used and corresponding files
# - create a directory per experiments and symlinks to corresponding files
# - create a compressed archive per dataset (following symlinks)

  # browser()
  fs::dir_create(export_dir, recurse = TRUE)
  
  # delete existing symlinks and directories to avoid creating undesired archive files
  if (length( fs::dir_ls(export_dir, type = c("directory", "symlink")) )) {
    if (!.force)
      stop('existing directories and/or symlinks in the export directory: please delete them or run with .force=TRUE (note that this option will delete them but force the compression of all files).')
  } else {
    system(sprintf('find %s -maxdepth 1 -type l -exec rm {} +', export_dir))
    system(sprintf('find %s -maxdepth 1 -mindepth 1 -type d -exec rm -rf {} +', export_dir))
  }


  myfiles %>% ungroup() %>% 
    # filter(!str_detect(condition, "stdIllum")) %>% 
    filter(str_detect(path, "data_thomas") | str_detect(path, "data_theo")) %>% 
    mutate(data_path = ifelse(str_detect(data_path, "^data"), paste0("./", data_path), data_path)) %>% 
    extract(data_path, "date_path", "(\\./data_[a-z]+/20\\d{6}/).*", remove = FALSE) %>% 
    select(date_path) %>% 
    distinct() %>%
    mutate(path = map(date_path, ~{#browser();
      fs::dir_ls(.) %>% as.character()
    })) %>%
    unnest(path) %>%
    filter(!str_detect(path, 'Pos\\d'), !str_detect(path, 'preproc'), !str_detect(path, 'output'), 
           !str_detect(path, 'curated'), !str_detect(path, '\\.properties')) %>%
    
    filter(
      !str_detect(path, "/data_theo/20210122/20210122_glu_lac_1$"),
      !str_detect(path, "20210122_short_test_script.sh"),
      !str_detect(path, "/data_theo/20210305/20210305_glu_lac_1$"),
      !str_detect(path, "20181008/old.bsh"),
      
    ) %>% 
    # arrange(path) %>%
    # pull(path)
    
    extract(date_path, "date", "\\./data_[a-z]+/(20\\d{6})/.*", remove = FALSE) %>% 
    mutate(filename = fs::path_file(path),
           epath = fs::path(export_dir, date, filename),
           is_dir = fs::is_dir(path)) %>%
    # print(n=Inf)
    # copy files and dirs
    # with( pwalk(list(path, epath), ~{ #browser();
    #   fs::dir_create(fs::path_dir(..2), recurse = TRUE)
    #   if (fs::is_file(..1)) fs::file_copy(..1, ..2)
    #   if (fs::is_dir(..1)) fs::dir_copy(..1, ..2)
    # }) ) %>%
    # create symlinks
    with( pwalk(list(path, epath), ~{ #browser();
      fs::dir_create(fs::path_dir(..2), recurse = TRUE)
      fs::link_create(fs::path_real(..1), ..2)
    }) ) %>%
    identity()
  
  # create files list before compressing
  # use bash command as fs doesnt allow to follow symlinks when traversing directories
  # # deprecated
  # system2("tree", c("-l", export_dir), stdout = fs::path(export_dir, paste0(fs::path_file(export_dir), '_fileList.txt')))
  # system2("find", c(export_dir, "-maxdepth 1 -mindepth 1 -type d -exec tree {} +"),
  #         stdout = fs::path(export_dir, paste0(fs::path_file(export_dir), '_fileList.txt')))
  # same but alphabetically sorted
  system2("find", c("-L", export_dir, "-mindepth 1 -type d -print0 | sort -z | xargs -r0 tree -l"),
          stdout = fs::path(export_dir , paste0(fs::path_file(export_dir), '_fileList.txt')))
  message("Done listing files...")

  # compress experiment by experiment (and resolve symlinks)
  fs::dir_walk(export_dir, function(.p) {# browser();
    if (fs::is_dir(.p) | (fs::is_link(.p))) {
      # tar(paste0(.p, '.tar.gz'), .p, compression = 'gzip', extra_flags = '--dereference')
      # fs::dir_delete(.p)
      # message("Done compressing ", .p)

      .f <- fs::path_file(.p)
      if (fs::file_exists(paste0(.p, '.tar.gz')) && !.force) {
        message ('skipping the archiving of ', .f)
        fs::dir_delete(.p)
      } else {
      sprintf("sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=tar_raw%s
#SBATCH --mem=32G
#SBATCH --time=1-0:00:00
#SBATCH --qos=1day
#SBATCH -o slogs/$JOB_NAME.o$JOB_ID
#SBATCH -e slogs/$JOB_NAME.e$JOB_ID
cd %s
tar --dereference -czf %s.tar.gz %s
rm -rf %s
EOF", .f, export_dir, .f, .f, .f) %>%
        system()
message("Job submitted for archiving ", .p)
    } }
  })

})()


(function(export_dir = "/scicore/home/nimwegen/GROUP/MM_Data/Thomas/_GCSlacArticle/Julou_2022_GCS_GL_Images",
          preproc_dir= "../../*preprocess*", .force=FALSE) {
  # this anonymous function exports preprocessed image data (each series of each experiment in separate compressed files) after the following scheme:
  # - traverse R data to list all datasets used
  # - create symlink to each dataset (renamed in a systematic manner)
  # - create a compressed archive per dataset (following symlinks)
  
  #    case of too many tif files: move them to another directory "toCurate" (already sorted by Theo)
  #    case of csv only (have been copied from the `preprocess` directory - write a script that find all export files 
  # recursively in a directory, and copy the parent subdirectory to the curated subfolder)

  # browser()
  fs::dir_create(export_dir, recurse = TRUE)
  
  # delete existing symlinks and directories to avoid creating undesired archive files
  if (length( fs::dir_ls(export_dir, type = c("directory", "symlink")) )) {
    if (!.force)
      stop('existing directories and/or symlinks in the export directory: please delete them or run with .force=TRUE (note that this option will delete them but force the compression of all files).')
  } else {
    system(sprintf('find %s -maxdepth 1 -type l -exec rm {} +', export_dir))
    system(sprintf('find %s -maxdepth 1 -mindepth 1 -type d -exec rm -rf {} +', export_dir))
  }
  
  
  # Thomas' data
  
  # start from myfiles rather than listing ./preproc in order to export only used data
  myfiles %>% ungroup() %>% 
    filter(str_detect(path, "data_thomas")) %>% 
    # filter(!str_detect(condition, "stdIllum")) %>% 
    mutate(data_path = ifelse(str_detect(data_path, "^data"), paste0("./", data_path), data_path)) %>% 
    select(path=data_path) %>% distinct() %>% 

    # copy files and dirs
    mutate(filename=fs::path_file(path)) %>% 
    extract(filename, 'date', '^(\\d{8})_', remove = TRUE) %>% 
    mutate(epath = fs::path(export_dir, date)) %>% 
    # View
    with( pwalk(list(path, epath), ~{ # browser();
      fs::dir_create(fs::path_dir(..2), recurse = TRUE)
      fs::link_create(fs::path_real(..1), ..2)
    }) ) %>%
    identity()
  
  
  # Theo's data - case: too many tif files
  
  # for d in series_dir:
  #   mkdir d in output_dir
  #   for f in tree(d)
  #     ln -s f output_dir/d
  #     for f2 in basenmae with tif but no csv
  #       rm f2
  
  myconditions %>% 
    map(~.$paths) %>% 
    unlist() %>% 
    tibble(path=.) %>% 
    filter(str_detect(path, "data_theo")) %>% 
    filter(str_detect(path, "20210122") | str_detect(path, "20210305")) %>% 

    # copy files and dirs
    mutate(filename=fs::path_file(path)) %>% 
    extract(filename, 'date', '^(\\d{8})_', remove = FALSE) %>% 
    extract(filename, 'series', 'S(\\d)/?$', remove = TRUE) %>% 
    mutate(ename = ifelse(is.na(series), date, paste(date, series, sep="_S")),
        epath = fs::path(export_dir, ename)) %>% 
    # print(n=Inf)
    with( pwalk(list(path, epath), ~{ # browser();
      fs::dir_create(fs::path_dir(..2), recurse = TRUE)
      fs::link_create(fs::path_real(..1), ..2)
    }) ) %>%
    identity()
  
  # delete all tif files for which no csv file is found
  vngMoM::find.files(export_dir, .follow_symlinks = T) %>% 
    tibble(path=.) %>% 
    extract(path, "basename", ".*/(?:ExportedCellStats_)?(20\\d{6}_.*_MMStack_Pos\\d{1,2}_preproc_GL\\d{2}).*", remove=FALSE) %>% 
    (function(.df) {
      anti_join(
        filter(.df, str_detect(path, ".tif$")),
        filter(.df, str_detect(path, ".csv$")),
        by="basename")
    }) %>% 
    pull(path) %>% 
    fs::file_delete()

  
  # Theo's data - case: csv only

  # for d in series_dir:
  #   mkdir d in output_dir
  #   for f in d:
  #     p = find f in preproc[d]
  #     ln -s dirname(p) output_dir/d
  
# tmp <-
  myconditions %>% 
    map(~.$paths) %>% 
    unlist() %>% 
    tibble(data_dir=.) %>% 
    filter(str_detect(data_dir, "data_theo")) %>% 
    filter(str_detect(data_dir, "20210504|20210506|20210513|20210708")) %>%
    mutate(path = map(data_dir, ~{#browser();
      fs::dir_ls(.) %>% as.character()
    })) %>%
    unnest(path) %>% 
    filter(!str_detect(path, "ExportedTracks")) %>% 
    extract(path, c("date", "pos", "gl"), "ExportedCellStats_(20\\d{6}).*Pos(\\d+)_GL(\\d+)\\.tiff?\\.csv", remove = FALSE) %>%
    extract(path, "gl_dir", "(Pos\\d+_GL\\d+)\\.tiff?\\.csv", remove = FALSE) %>% 
    # identity()
    # tmp2 <- tmp %>%
    left_join(
      # Warning: there can be several preprocess directories... look only once per position for the correct one
      .,
      group_by(., date, pos) %>% 
        slice(1L) %>% 
        mutate(source_dir=pmap_chr(
          list(data_dir, path), 
          ~ vngMoM::find.files(fs::path(..1, preproc_dir), .name=basename(..2)) %>% 
            dirname() %>% fs::path("..") %>% fs::path_norm()
        )) %>% 
        select(date, pos, source_dir),
      by=c("date", "pos")) %>% 
    mutate(source_path=fs::path(source_dir, paste0("Pos", pos, "_GL", gl)),
           epath = fs::path(export_dir, basename(data_dir), paste0("Pos", pos, "_GL", gl))) %>% 
    # print(n=Inf)
    with( pwalk(list(source_path, epath), ~{ # browser();
      fs::dir_create(fs::path_dir(..2), recurse = TRUE)
      fs::link_create(fs::path_real(..1), ..2)
    }) ) %>%
    identity()
  
  
  # create files list before compressing
  # use bash command as fs doesnt allow to follow symlinks when traversing directories
  # # deprecated
  # system2("tree", c("-l", export_dir), stdout = fs::path(export_dir, paste0(fs::path_file(export_dir), '_fileList.txt')))
  # system2("find", c(export_dir, "-maxdepth 1 -type l -print0 | sort -z | xargs -r0 tree -l"), 
  #         stdout = fs::path(export_dir, paste0(fs::path_file(export_dir), '_fileList.txt')))
  system(paste0(
    "(find ", export_dir, " -maxdepth 1 -type l -print0; find ", 
    export_dir, " -mindepth 1 -type d -print0) | sort -z | xargs -r0 tree -l > ",
    fs::path(export_dir , paste0(fs::path_file(export_dir), '_fileList.txt')) ))
  message("Done listing files...")
  
  # compress experiment by experiment (and resolve symlinks)
  fs::dir_walk(export_dir, function(.p) {# browser();
    if (fs::is_dir(.p) | (fs::is_link(.p))) {
      # tar(paste0(.p, '.tar.gz'), .p, compression = 'gzip', extra_flags = '--dereference')
      # fs::dir_delete(.p)
      # message("Done compressing ", .p)
      .f <- fs::path_file(.p)
      if (fs::file_exists(paste0(.p, '.tar.gz')) && !.force) {
        message ('skipping the archiving of ', .f)
        try(fs::link_delete(.p))
        try(fs::dir_delete(.p))
      } else {
        sprintf("sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=tar_GL%s
#SBATCH --mem=32G
#SBATCH --time=1-0:00:00
#SBATCH --qos=1day
#SBATCH -o slogs/$JOB_NAME.o$JOB_ID
#SBATCH -e slogs/$JOB_NAME.e$JOB_ID
cd %s
tar --dereference -czf %s.tar.gz %s
rm -rf %s
EOF", .f, export_dir, .f, .f, .f) %>%
        system()
      # message("Job submitted for archiving ", .p)
    } }
  })

})()


                                                                                                                                                                                                            (function(export_dir = "/scicore/home/nimwegen/GROUP/MM_Data/Thomas/_GCSlacArticle/Julou_2022_GCS_GL_Preproc") {
  # this anonymous function exports parsed output (tabular data) after the following scheme:
  # - traverse R data to list all datasets used
  # - create symlink to each dataset (renamed in a systematic manner)
  # - create a compressed archive (following symlinks)
  
  # browser()
  fs::dir_create(export_dir, recurse = TRUE)
  
  # delete existing symlinks and directories to avoid creating undesired archive files
  if (length( fs::dir_ls(export_dir, type = c("directory", "symlink")) )) {
    if (!.force)
      stop('existing directories and/or symlinks in the export directory: please delete them or run with .force=TRUE (note that this option will delete them but force the compression of all files).')
  } else {
    system(sprintf('find %s -maxdepth 1 -type l -exec rm {} +', export_dir))
    system(sprintf('find %s -maxdepth 1 -mindepth 1 -type d -exec rm -rf {} +', export_dir))
  }
  
  # start from myfiles rather than listing ./preproc in order to export only used data and to split date by series
  myfiles %>% ungroup() %>% 
    select(cpath=path) %>% 
    mutate(cdir = fs::path_dir(cpath), cname=fs::path_file(cpath)) %>% 
    filter(str_detect(cpath, "data_thomas") | str_detect(cpath, "data_theo")) %>% 
    extract(cpath, 'date', 'data_th[a-z]+/(\\d{8})/', remove = FALSE) %>% 
    extract(cdir, 'series', 'S(\\d)/?$', remove = FALSE) %>% 
    mutate(edir = ifelse(is.na(series), date, paste(date, series, sep="_S")),
           path = data2preproc(cpath),
           epath = fs::path(export_dir, edir, data2preproc_file(cname)) ) %>%
  # select(path, epath) %>% print(n=Inf)
  # myconditions %>%
  #   map(~.$paths) %>% 
  #   unlist() %>% 
  #   tibble(cpath=.) %>% 
  #   filter(str_detect(cpath, "data_thomas") | str_detect(cpath, "data_theo")) %>% 
  #   extract(cpath, 'date', 'data_th[a-z]+/(\\d{8})/', remove = FALSE) %>% 
  #   extract(cpath, 'series', 'S(\\d)/?$', remove = FALSE) %>% 
  #   mutate(ename = ifelse(is.na(series), date, paste(date, series, sep="_S")),
  #          path = fs::path("preproc", ename),
  #          epath = fs::path(export_dir, ename) ) %>% 
  # create symlinks
  with( pwalk(list(path, epath), ~ {
    fs::dir_create(fs::path_dir(..2), recurse = TRUE)
    fs::link_create(fs::path_real(..1), ..2, ) 
  }) ) %>%
    identity()
  
  # create files list before compressing
  system2("find", c("-L", export_dir, "-mindepth 1 -type d -print0 | sort -z | xargs -r0 tree"),
          stdout = fs::path(export_dir, '..' , paste0(fs::path_file(export_dir), '_fileList.txt')))
  message("Done listing files...")
  
  # system2("tree", c("-l", export_dir), stdout = fs::path(export_dir, '..', paste0(fs::path_file(export_dir), '_fileList.txt')))
  
  # compress all together (and resolve symlinks)
  sprintf("cd %s; tar --dereference -czf %s.tar.gz %s", 
          fs::path_dir(export_dir), fs::path_file(export_dir), fs::path_file(export_dir)) %>% 
    system()

  fs::dir_delete(export_dir)
})()



myframes %>% ungroup() %>% 
  select(path, condition, date) %>% 
  distinct() %>% 
  filter(str_detect(path, "data_thomas") | str_detect(path, "data_theo")) %>% 
  extract(path, 'series', '20\\d{6}_S(\\d)/', remove = FALSE) %>% 
  mutate(edir = ifelse(is.na(series), date, paste(date, series, sep="_S"))) %>% 
  select(dir=edir, condition) %>% 
  distinct() %>% 
  filter(!str_detect(condition, "stdIllum")) %>% 
  arrange(dir) %>% 
  write_excel_csv("/scicore/home/nimwegen/GROUP/MM_Data/Thomas/_GCSlacArticle/Julou_2022_GCS_series.csv")
  
         
         
# # IDR metadata #######
# (function() {
#   # this anonymous function exports IDR metadata after the following scheme:
#   # - traverse R data to list all datasets used and corresponding files
#   # - create a directory per experiments and symlinks to corresponding files
#   # - create a compressed archive per dataset (following symlinks)
#   
#   bind_rows(
#     myconditions %>% 
#       map(~{tibble(condition=.$condition, cpath=.$paths)}) %>% 
#       bind_rows() %>% 
#       filter(str_detect(cpath, "data_thomas")) %>% 
#       extract(cpath, 'date', 'data_thomas/(\\d{8})/', remove = FALSE) %>% 
#       mutate(path = map(cpath, ~{#browser();
#         fs::dir_ls(fs::path(., '..')) %>% as.character()
#       })) %>% 
#       unnest(path) %>% 
#       mutate(path=fs::path_norm(path)) %>% 
#       
#       filter(!str_detect(path, 'Pos\\d'), !str_detect(path, 'curated')) %>%
#       filter(!str_detect(path, '20151127_flatfield'),
#              !str_detect(path, '20190605_ASC662'),
#              !str_detect(path, '20170108_gluLac_lac_snapAfter'),
#              !path %in% c('data_thomas/20180514/20180515_gfpFlatfield_pos',
#                           'data_thomas/20180606/20180606_gluLac_lac_switch16h_1')) %>% 
#       filter(fs::is_dir(path)) %>% 
#       mutate(path=map_chr(path, ~fs::dir_ls(., glob="*.ome.tif")[1] ),
#              file=fs::path_file(as.character(path))
#       ) %>% 
#       # with(interaction(date, file)) %>% unique %>% length
#       mutate(type=ifelse(str_detect(path, fixed('flatfield', ignore_case=TRUE)), 'flatfield', 'raw')) %>% 
#       select(condition, dataset=date, file, path, type) %>% 
#       mutate(path=str_replace(path, 'data_thomas', 'Julou_2020_lacInduction_RawImages'),
#              preproc='') %>% 
#       identity(),
#     
#     myconditions %>% 
#       map(~{tibble(condition=.$condition, cpath=.$paths)}) %>% 
#       bind_rows() %>% 
#       filter(str_detect(cpath, "data_thomas")) %>% 
#       mutate(filename=fs::path_file(cpath)) %>% 
#       extract(filename, 'date', '^(\\d{8})_', remove = TRUE) %>%
#       
#       # filter(date %in% c('20151204', '20180313')) %>%
#       
#       mutate(path=map(cpath, ~fs::dir_ls(., type = "directory", recurse=T) %>% as.character()) ) %>% 
#       unnest(path) %>% 
#       mutate(path=map_chr(path, ~fs::dir_ls(., glob="*.tif")[1]),
#              file=fs::path_file(as.character(path)) ) %>% 
#       mutate(type='derived') %>% 
#       extract(path, c("pos", "gl"), ".*\\d{8}_.*[Pp]os(\\d+).*_GL(\\d+).*", remove=FALSE) %>%
#       left_join(
#         myfiles %>% 
#           mutate(preproc=data2preproc(path)) %>% 
#           extract(path, c("date", "pos", "gl"), ".*(\\d{8})_.*[Pp]os(\\d+).*_GL(\\d+).*", remove=FALSE) %>%
#           ungroup() %>% select(date, pos, gl, preproc) %>% 
#           identity(),
#         by=c('date', 'pos', 'gl')
#       ) %>% 
#       select(condition, dataset=date, file, path, type, preproc) %>% 
#       mutate(path=str_replace(path, '\\./data_thomas', 'Julou_2020_lacInduction_GL_Images')) %>% 
#       mutate(preproc=str_replace(preproc, '\\./preproc', 'GL_Preproc')) %>% 
#       identity(),
#   ) %>% 
#     
#     left_join(readxl::read_xlsx('data_thomas/MM_Experiments_Table.xlsx') %>% 
#                 transmute(dataset=date, strain=str_replace_all(strain, '\r?\n', ' | '), 
#                           media = str_replace_all(media, '\r?\n', ' | '), 
#                           steps = str_replace_all(steps, '\r?\n', ' | '), interval_min, 
#                           flow_control = str_replace_all(flow_control, '\r?\n', ' | ')),
#               by='dataset') %>% 
#               
#     transmute(
#       'Source Name' = 'Julou_2020-lacInduction',  # an identifier for the sample
#       'Characteristics [Organism]' = 'Escherichia coli',  # the species of the sample e.g.Homo sapiens
#       'Term Source 1 REF' = 'NCBITaxon',
#       'Term Source 1 Accession' = '', # leave blank
#       'Characteristics [Cell Line]' = 'ASC662 (MG1655 lacZ-GFPmut2)', # the name of the cell line used e.g HeLa (if appropriate)  Other Characteristics columns can be added as necessary e.g. Characteristics [Organism Part]
#       'Term Source 2 REF' = 'EFO',
#       'Term Source 2 Accession' = '', # leave blank
#       'Protocol REF' = 'growth protocol',
#       'Protocol REF ' = 'treatment protocol', 
#       'Protocol REF  ' = 'image acquisition and feature extraction protocol',
#       'Protocol REF   ' = 'image analysis protocol', # data analysis protocol
#       'Assay Name' = condition,  # an name for the imaging assay. This column can be used to group several images together If several images (e.g. fields, or replicate images) have been taken of the same sample then repeat the assay name on each row corresponding to the imaging file.  It can also be used to group raw and processed images from the same assay together
#       'Experimental Condition [Media]' = media, # add an experimental condition if there is one e.g. 'targeted protein', or 'antibody'
#       'Experimental Condition [Steps]' = steps, # add an experimental condition if there is one e.g. 'targeted protein', or 'antibody'
#       'Experimental Condition [Flow Control]' = flow_control, # add an experimental condition if there is one e.g. 'targeted protein', or 'antibody'
#       'Comment [Preculture]' = strain, 
#       # 'Comment [Acquisition Interval]' = interval_min, 
#       'Comment [Gene Identifier]' = '', # enter an identifier for any associated genes e.g. ensembl identifier for targeted protein
#       'Comment [Gene Symbol]' = '', # enter gene symbol for any associated genes e.g. BRCA1
#       'Comment [Gene Annotation Comments]' = '', # any comments about the gene annotation e.g. what gene annotation build the gene identifiers come from
#       'Dataset Name' = dataset, # this column can be used to group images into datasets. These datasets will be used to group the images in the Image Data Resource. The name can be the same as the assay name.
#       'Image File' = file,  # the name of the image file
#       'Comment [Image File Path]' = path, # the path to the file
#       'Comment [Image File Comments]' = '', # any comments about image files (can say its missing if you no longer have it)
#       'Comment [Image File Type]' = type, # is it a raw image or derived?
#       'Channels' = "Phase Contrast, GFP [epifluorescence]",  # the names of the channels and what was labeled in each channel
#       'Processed Data File' = preproc, # name of the file with the results in
#     ) %>% 
#     
#   write_delim('share/Julou_2020-lacInduction-assays.txt', delim='\t')
# 
# })()
# 
# 
# 
# # PLoS (lags export) #######
# 
# mycells_switching %>% ungroup %>% 
#   filter_article_ds() %>% 
#   filter(! date %in% discarded_dates) %>% 
#   filter(!discard_arrested | condition=='switch_gly_lac') %>%
#   filter(switch_idx == 1 | (str_detect(condition, '^switch_[0-9]+h$') & switch_idx < 3) ) %>% 
#   left_join(
#     mycells_switching_memory %>% 
#       ungroup() %>% 
#       mutate(gfp_inherit=parent_gfp/2^divs_since_par) %>% 
#       select(ugen=ugen.x, gfp_inherit, parent_gfp, divs_since_par)
#   ) %>% 
#   left_join(select(mycells_switchingLow_memory, ugen, gfp_inherit, fluo_focus)) %>% 
#   select(condition, ugen, date, pos, gl, switch_idx, lag_gfp_200=lag_200, growth_lag, 
#          gfp_ini, gr_before=logl_time_slope_before, gfp_inherit, parent_gfp, n_divs_since_parent=divs_since_par, fluo_focus, time_birth, time_div) %>% 
#   write.csv("share/Julou_2020-lacInduction-lags.csv", row.names=FALSE)
#   # identity()
# 

