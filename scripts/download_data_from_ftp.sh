# make /data folder and subfolders
mkdir -p /ibangs
mkdir -p /ibangs/data
mkdir -p /ibangs/tutorial

# Copy the data directory.
wget --directory-prefix=/ibangs/data ftp://ftp.jax.org/dgatti/IBANGS2016/data/DO_Sanger_SDPs.txt.bgz
wget --directory-prefix=/ibangs/data ftp://ftp.jax.org/dgatti/IBANGS2016/data/DO_Sanger_SDPs.txt.bgz.tbi
wget --directory-prefix=/ibangs/data ftp://ftp.jax.org/dgatti/IBANGS2016/data/assoc_perms.rds
wget --directory-prefix=/ibangs/data ftp://ftp.jax.org/dgatti/IBANGS2016/data/expr_sig_qtl.rds
wget --directory-prefix=/ibangs/data ftp://ftp.jax.org/dgatti/IBANGS2016/data/haploprobs.rds
wget --directory-prefix=/ibangs/data ftp://ftp.jax.org/dgatti/IBANGS2016/data/ibangs_expr.rds
wget --directory-prefix=/ibangs/data ftp://ftp.jax.org/dgatti/IBANGS2016/data/ibangs_haploprobs.rds
wget --directory-prefix=/ibangs/data ftp://ftp.jax.org/dgatti/IBANGS2016/data/ibangs_phenotypes.rds
wget --directory-prefix=/ibangs/data ftp://ftp.jax.org/dgatti/IBANGS2016/data/linkage_perms.rds
wget --directory-prefix=/ibangs/data ftp://ftp.jax.org/dgatti/IBANGS2016/data/phenotypes.rds

# Copy the tutorial directory.
wget --directory-prefix=/ibangs/tutorial ftp://ftp.jax.org/dgatti/IBANGS2016/markdown/DO.impute.founders.sm.png
wget --directory-prefix=/ibangs/tutorial ftp://ftp.jax.org/dgatti/IBANGS2016/markdown/DOQTL_workshop_IBANGS2016.Rmd
wget --directory-prefix=/ibangs/tutorial ftp://ftp.jax.org/dgatti/IBANGS2016/markdown/DOQTL_workshop_IBANGS2016.html

# set privilages - everybody can do everything
chmod --recursive 777 /ibangs
