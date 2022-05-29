library(data.table)

read_path = '/home/panwei/shared/UKBiobankIndiv/ukb49020.tab'
write_path = '~/deepRIV/UKB/data/'
# covariates
cov_col = c('f.eid', 'f.21003.0.0', 'f.22001.0.0',
            'f.22009.0.1','f.22009.0.2','f.22009.0.3','f.22009.0.4','f.22009.0.5',
            'f.22009.0.6','f.22009.0.7','f.22009.0.8','f.22009.0.9','f.22009.0.10')
cov_names = c('f.eid','age','sex',
              'pc1','pc2','pc3','pc4','pc5',
              'pc6','pc7','pc8','pc9','pc10')

covariates = fread(read_path, select = cov_col)
colnames(covariates) = cov_names
fwrite(covariates, paste0(write_path,'covariates.csv'))

# exposures
exp_col = c('f.eid', 'f.23104.0.0', 'f.23099.0.0', 'f.48.0.0', 'f.49.0.0')
exp_names = c('f.eid', 'bmi', 'bfp', 'waist-circum', 'hip-circum')
exposures = fread(read_path, select = exp_col)
colnames(exposures) = exp_names
fwrite(exposures, paste0(write_path,'exposures.csv'))

# outcomes
out_col = c('f.eid', 'f.102.0.0', 'f.4080.0.0', 'f.4079.0.0')
out_names = c('f.eid', 'pulse-rate', 'sbp', 'dbp')
outcomes = fread(read_path, select = out_col)
colnames(outcomes) = out_names
fwrite(outcomes, paste0(write_path,'outcomes.csv'))
