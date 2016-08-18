## SNMM_Start_Job_By_ID(filename_job_specifications,job_ID)
SNMM_Stan_Start_Job_By_ID <- function(filename_job_specifications,job_ID)
{
  require("rstan")
  job_specs <- read.csv(file=filename_job_specifications,header=TRUE,sep=";")
  print(job_specs)
  fn_datafile <- toString(job_specs$filename[job_ID])
  data_snmm_thinned <- SNMM_Load_Data(fn_datafile,as.numeric(job_specs$data_thinning_factor[job_ID]))
  
  ## plot a histogram of the data (just for temporary test purposes:)
  hist(data_snmm_thinned,breaks=100,main="Histogram of deflection data",xlab="deflection [nm]",ylab="frequency")
  ##
  
  nstates <- as.numeric(job_specs$num_states[job_ID])
  num_samples <- as.numeric(job_specs$num_total_samples[job_ID])
  sample_thinning_factor <- as.numeric(job_specs$sample_thinning_factor[job_ID])
  num_chains <- as.numeric(job_specs$num_chains[job_ID])
  locations0 <- as.numeric(job_specs[job_ID,7:(7+nstates-1)])
  scales0 <- as.numeric(job_specs[job_ID,(7+nstates):(7+2*nstates-1)])
  shapes0 <- as.numeric(job_specs[job_ID,(7+2*nstates):(7+3*nstates-1)])
  pis0 <- as.numeric(job_specs[job_ID,(7+3*nstates):(7+4*nstates-1)])
  alphas <- array(1,dim=c(nstates))
  
  print(locations0)
  print(scales0)
  print(shapes0)
  print(pis0)
  
##  snorm.analysis <- smsn.mix(y=data_snmm_thinned,mu=locations0,sigma2=scales0,shape=shapes0,pii=pis0,g=2,nu=3,get.init=FALSE,group=TRUE,family="Skew.normal")
  
  
  num_points <- length(data_snmm_thinned)
  data_stan_snmm <- list(N=num_points,M=nstates,y=data_snmm_thinned,alphas=alphas,mus0=locations0,std0=2);
  init_stan_snmm <- list(list(phis=pis0,locations=locations0,scales=scales0,shapes=shapes0))
  stan_fit_snmm <- stan(file="snmm_general.stan",data=data_stan_snmm,init=init_stan_snmm,iter=num_samples,chain=num_chains,thin=sample_thinning_factor)
  
  mat_stan_snmm_samples <- as.matrix(stan_fit_snmm)
  ## obtain the "case-string":
  splitted <- strsplit(x=fn_datafile,split=".",fixed=TRUE)
  casename <- splitted[[1]][1]
  
  filename_snmm_samples <- paste0(casename,"_snmm_samples.csv")
  write.table(x=mat_stan_snmm_samples,file=filename_snmm_samples,sep=";")
  SNMM_Stan_Evaluate_Samples(stan_fit_snmm,data_snmm_thinned,num_states=nstates,casename=casename)
  return(stan_fit_snmm)
}

## SNMM_Stan_Evaluate_Samples(stan_fit_snmm,num_states=2,casename="test_case")
## input: takes a stan-fit object "stan_fit_snmm", the original data "data", the
## number of states "num_states" (by default 2) and a casename "casename" (by default: "test_case")
## as an input. 
## output: plots a histogram together with the model (mean-posterior parameters)
## and writes the values of the mean-posteriors to a file
SNMM_Stan_Evaluate_Samples <- function(stan_fit_snmm,data,num_states=2,casename="test_case")
{
  filename_snmm_plot <- paste0(casename,"_snmm_plot.pdf")
  pdf(file=filename_snmm_plot,width=14,height=11,useDingbats=FALSE)
  SNMM_Stan_Plot_Skewed_Normals_From_Mean_Posterior(data,stan_fit_snmm,num_states,case_string=casename)
  dev.off()
  
  phi_means <- array(0,dim=c(num_states))
  location_means <- array(0,dim=c(num_states))
  scale_means <- array(0,dim=c(num_states))
  shape_means <- array(0,dim=c(num_states))
  
  mat <- as.matrix(stan_fit_snmm)
  for(i in 1:num_states)
  {
    phi_means[i] <- mean(mat[,i]);
    location_means[i] <- mean(mat[,i+num_states]);
    scale_means[i] <- mean(mat[,i+2*num_states]);
    shape_means[i] <- mean(mat[,i+3*num_states]);
  }
  
  df_posterior_means <- data.frame(phis=phi_means,locs = location_means,scales=scale_means,shapes=shape_means)
  
  filename_posterior_mean <- paste0(casename,"_snmm_post_mean.csv")
  write.table(x=df_posterior_means,file=filename_posterior_mean,sep=";")
}

## SNMM_Stan_Plot_Skewed_Normals_From_Mean_Posterior(data_orig,fit_stan_snmm,
## num_states=2,lab_x="deflection [nm]",case_string)
## input: the original data "data_orig" and a stan-fit object
## "fit_stan_snmm" (skewed normal mixture model-fit), as well as the
## number of states "num_states", the label of the x-axis and 
## a case_string "case_string"
## output: plots a histogram of the original data, as well as
## the mixture model (components individually as well as combined)
SNMM_Stan_Plot_Skewed_Normals_From_Mean_Posterior <- function(data_orig,fit_stan_snmm,num_states=2,lab_x="deflection [nm]",case_string)
{
  mat <- as.matrix(fit_stan_snmm);
  phis_means <- matrix(nrow=num_states,ncol=1);
  location_means <- matrix(nrow=num_states,ncol=1);
  scale_means <- matrix(nrow=num_states,ncol=1);
  shape_means <- matrix(nrow=num_states,ncol=1);
  
  for(i in 1:num_states)
  {
    phis_means[i] <- mean(mat[,i]);
    location_means[i] <- mean(mat[,i+num_states]);
    scale_means[i] <- mean(mat[,i+2*num_states]);
    shape_means[i] <- mean(mat[,i+3*num_states]);
  }
  SNMM_Stan_Plot_Skewed_Normals_Mixture_Model(data_orig,phis_means,location_means,scale_means,shape_means,lab_x,case_string);
}

## PlotSkewedNormalsMixtureModel(data_orig,phis,mus,sigmas,alphas)
## input: the original (either force or extension) data "data_orig",
## the state probabilities "phis", the mean values "mus", the standard-deviations
## "sigmas" and the scales "alphas"
SNMM_Stan_Plot_Skewed_Normals_Mixture_Model <- function(data_orig,phis,locations,scales,shapes,lab_x="force [pN]",case_string)
{
  require("sn")
  num_states <- length(phis)
  minval <- min(data_orig)
  maxval <- max(data_orig)
  delta <- (maxval - minval)/200
  minval <- minval-10*delta
  maxval <- maxval+10*delta
  xs <- seq(from=minval,to=maxval,by=delta)
  ys_temp <- matrix(nrow=length(xs),ncol=1)
  ys <- matrix(0,nrow=length(xs),ncol=1)
  hist(data_orig,breaks=200,freq=FALSE,xlab=lab_x,xlim=c(minval,maxval))
  color_table <- c("blue","red","green","yellow","cyan")
  for(j in 1:num_states)
  {
    for(i in 1:length(xs))
    {
      ys_temp[i] <- phis[j]*dsn(xs[i],xi=locations[j],omega=scales[j],alpha=shapes[j])
    }
    points(xs,ys_temp,type="l",col=color_table[j],lwd=2)
    ys <- ys + ys_temp
  }
  points(xs,ys,type="l",col="black",lwd=2)
}

## SNMM_Load_Data(filename,thinning_factor)
## input: 
##
SNMM_Load_Data <- function(filename,thinning_factor)
{
  data_snmm <- read.table(file=filename,header=TRUE)
  data_snmm_thinned <- data_snmm[seq(from=1,to=length(data_snmm[,1]),by=thinning_factor),1]
  return(data_snmm_thinned)
}