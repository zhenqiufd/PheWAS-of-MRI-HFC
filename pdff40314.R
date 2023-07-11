###PheWAS of hepatic fat content###
###start time: 2023-04-26       ###
###author: Zhenqiu Liu          ###
###---------Linux for PDFF GWAS-----------------------------------------------------------
library(tidyverse)
library(data.table)
library(doParallel)
library(bigsnpr)
library(bigreadr)
setwd('/disk/disk1/UKB')
###GWAS FOR PDFF (n=40314)
pdff <- read.csv('/share/home/zhenqiu/PDFF_40314.csv')
#install.packages('RNOmni')
library(RNOmni)
pdff2 <- pdff %>% select(1,3) %>%
  rename(eid = eid.63726) %>%
  filter(!is.na(pdff)) %>%
  mutate(FRPDFF = RankNorm(pdff, k = 3/8))

registerDoParallel(cl <- makeCluster(30))
system.time({
  list_snp_id <- foreach(chr = 1:22) %dopar% {
    # cat("Processing chromosome", chr, "..\n")
    mfi <- paste0("/disk/disk1/UKB/Gene/Imputation/ukb_mfi_chr", chr, "_v3.txt")
    infos_chr <- bigreadr::fread2(mfi, showProgress = FALSE)
    infos_chr_sub <- dplyr::filter(infos_chr,V6 > 0.01, V8 > 0.5)
    #V8 denotes the imputation quality; V6 MAF
    with(infos_chr_sub, paste(chr, V3, V4, V5, sep = "_"))
  }
}) 
stopCluster(cl)

sample <- fread2("/disk/disk1/UKB/Gene/Imputation/sampleID/ukb63726_imp_chr1_v3_s487296.sample")
str(sample)

(N <- readBin("/disk/disk1/UKB/Gene/Imputation/ukb_imp_chr1_v3.bgen", what = 1L, size = 4, n = 4)[4])
sample <- sample[-1, ]
stopifnot(nrow(sample) == N)

csv <- "/disk/disk1/UKB/Data/ukb_data.csv"
df0 <- fread2(csv, select = c("f.eid","f.31.0.0","f.34.0.0","f.53.0.0",
                              "f.22001.0.0","f.21000.0.0","f.21001.0.0",
                              "f.22006.0.0","f.22027.0.0"))
colnames(df0) <- c("f.eid","reportedSex","YearOfBirth","attendingDate",
                   "geneticSex","EthnicBackground",'BMI',
                   "is_caucasian","isOutliers")
df1 <- df0 %>% filter(is.na(isOutliers),
                      reportedSex==geneticSex,
                      is_caucasian==1)

ind.indiv <- match(df1$f.eid, sample$ID_2)
sub <- which(!is.na(ind.indiv))
sub_eid <- df1$f.eid[sub]
df2 <- df1 %>% filter(f.eid %in% sub_eid)

# Relatedness
rel <- fread2("/disk/disk1/UKB/Gene/Relatedness/ukb63726_rel_s488264.dat")
sum(ind.rm <- df2$f.eid %in% rel[rel$Kinship > 0.08, "ID2"])
df3 <- df2[!ind.rm, ]

df4 <- df3 %>% left_join(pdff2,by=c('f.eid'='eid')) 
df5 <- df4 %>% filter(!is.na(FRPDFF))
#n=31377
#Genotypes
ind.indiv2 <- match(df5$f.eid, sample$ID_2)
sub2 <- which(!is.na(ind.indiv2))

system.time(
  rds <- bigsnpr::snp_readBGEN(
    bgenfiles = glue::glue("/disk/disk1/UKB/Gene/Imputation/ukb_imp_chr{chr}_v3.bgen", 
                           chr = 1:22),
    list_snp_id = list_snp_id,
    backingfile = "/disk/disk1/workspace/liuzhenqiu/UKBB_full_V8_50_pdff",
    ind_row = ind.indiv2[sub2],
    #bgi_dir = "Gene/Imputation",
    ncores = 50
  )
) # 26.3 mins

simu <- snp_attach("/disk/disk1/workspace/liuzhenqiu/UKBB_full_V8_50_pdff.rds")
G <- simu$genotypes
length(G)

pcs <- fread2(csv,select = c('f.eid',paste0('f.22009.0.',1:10)))
pcs <- pcs %>% rename(eid=f.eid, sex = 'f.31.0.0',
                      yearofbirth='f.34.0.0',
                      pc1 = f.22009.0.1,
                      pc2 = f.22009.0.2,pc3 = f.22009.0.3,
                      pc4 = f.22009.0.4,pc5 = f.22009.0.5,
                      pc6 = f.22009.0.6,pc7 = f.22009.0.7,
                      pc8 = f.22009.0.8,pc9 = f.22009.0.9,
                      pc10 = f.22009.0.10)
df6 <- df5 %>% left_join(pcs, by = c('f.eid'='eid'))
#BMI <- ifelse(is.na(df6$BMI),median(df6$BMI,na.rm=T),df6$BMI)
age2 <- (df6$age)^2
COVAR <- cbind(df6$yearofbirth,df6$sex,age2,
               df6$pc1,df6$pc2,df6$pc3,df6$pc4,df6$pc5,
               df6$pc6,df6$pc7,df6$pc8,df6$pc9,df6$pc10)
#COVAR <- apply(COVAR,2,function(x) ifelse(is.na(x),median(x,na.rm=T),x))
y <- df6$FRPDFF
system.time(
  gwas0 <- big_univLinReg(G, y, covar.train = COVAR, ncores = 1)
) 

ggsave(filename = 'snp_qq_pdff.jpeg',
       plot = snp_qq(gwas0),
       dpi = 300,
       height = 15,
       width = 15,
       units = 'cm',
       path = '/disk/disk1/workspace/liuzhenqiu'
)
gwas0$p.value <- predict(gwas0, log10 = FALSE)
data.table::fwrite(gwas0,"/disk/disk1/workspace/liuzhenqiu/gwas0_PDFF_mhtest.csv")
#gwas0$p.value_log <- log(gwas0$p.value)
#fivenum(gwas0$p.value_log)
sum(gwas0$p.value< 5e-8)
##n=939
snp_ids <- unlist(list_snp_id)
gwas0 <- gwas0 %>% mutate(snp_id = snp_ids)
gwas_sig <- gwas0 %>% filter(p.value< 5e-8) %>% arrange(p.value)

rsid <- data.table::fread("/disk/disk1/workspace/liuzhenqiu/snp_rsid_hg19.txt.gz")
#sigSNP <- fread('gwas0_PDFF_sig0318.csv')
gwas1 <- gwas_sig %>% mutate(
  chr = str_extract(snp_id,'^[0-9]{1,}'),
  start = str_extract(snp_id,'[0-9]{5,}'),
  `chromosome:start` = paste0(chr,":",start)
)

gwas1 <- gwas1 %>% left_join(rsid,by = 'chromosome:start')
gwas1 <- gwas1 %>% rename(beta = estim,
                          se = std.err,
                          pval = p.value,
                          pos = start,
                          snp = name)
saveRDS(gwas1,'/disk/disk1/workspace/liuzhenqiu/sig_snps_pdff.rds')

gwas0 <- gwas0 %>% mutate(
  chr = str_extract(snp_id,'^[0-9]{1,}'),
  start = str_extract(snp_id,'[0-9]{5,}'),
  `chromosome:start` = paste0(chr,":",start)
)
gwas0 <- gwas0 %>% left_join(rsid,by = 'chromosome:start')
gwas0 <- gwas0 %>% rename(beta = estim,
                          se = std.err,
                          pval = p.value,
                          pos = start,
                          snp = name)
saveRDS(gwas0,'/disk/disk1/workspace/liuzhenqiu/gwas0_PDFF_mhtest0426.rds')
library(ieugwasr)
x = ld_clump(tibble(rsid=gwas1$snp,pval=gwas1$pval),
             clump_r2=0.01,
             plink_bin='/share/home/zhenqiu/plink',
             bfile='/share/home/zhenqiu/EUR-LD/EUR')
independent_sig_snps <- gwas1 %>% filter(snp %in% x$rsid)
write.csv(independent_sig_snps,'/disk/disk1/workspace/liuzhenqiu/independent_sig_snps.csv',row.names = F)

####query the snps in linux
library(rbgen)
query_snps <- read.csv('/disk/disk1/workspace/liuzhenqiu/independent_sig_snps.csv')
sample2 <- sample
snps <- query_snps$snp
for(i in 1:length(snps)){
  snp = snps[i]
  chr = query_snps$chr[i]
  dirs = paste0('/disk/disk1/UKB/Gene/Imputation/ukb_imp_chr',chr,'_v3.bgen')
  index.filename = paste0('/disk/disk1/UKB/Gene/Imputation/ukb_imp_chr',chr,'_v3.bgen.bgi')
  snp_info = bgen.load(dirs,rsids = snp,index.filename=index.filename)
  rsid <- as_tibble(snp_info$data[,,1:3])
  try1 <- try({sample2 <- sample2 %>% bind_cols(rsid)},silent = T)
  if(inherits(try1,'try-error')) {
    print(paste0(i,' was skipped!!!'))
    next
  } else {
    sample2 <- sample2 %>% mutate(snp = ifelse(round(`g=0`) == 1,0,ifelse(round(`g=1`)==1,1,2)))
    sample2 <- sample2 %>% select(-c(`g=0`,`g=1`,`g=2`))
    names(sample2)[ncol(sample2)] <- snps[i]
    print(paste0(i,' was finished!!!'))
  }}
fwrite(sample2,'/disk/disk1/workspace/liuzhenqiu/snps_query_data.csv')

####----------WINDOWS for subsequent analyses------------------------------------------------------------
###---manhattan plot-------------------------------------
setwd('E:\\English Papers for writing\\PDFF-GWAS-MR')
gwasPDFF <- readRDS('gwas0_PDFF_mhtest0426.rds')
source('E:\\manhattanPlot.R')
pdff2 <- gwasPDFF %>% select(pval,chr,pos,snp) %>%
  rename(SNP=snp,CHR=chr,P=pval,BP=pos)
pdff2 <- pdff2 %>% mutate(
  CHR = as.numeric(CHR),
  BP = as.numeric(BP)
)
ggsave(
  filename = 'PDFF_manhattan.jpeg',
  myManhattan(pdff2,y.step=20,suggestiveline=F),
  dpi=300,
  width = 45,
  height = 30,
  units = 'cm'
)

###---phenotypic data processing--------------------------------------------------------
setwd('E:\\English Papers for writing\\PDFF-GWAS-MR')
library(tidyverse)
library(data.table)
library(lubridate)
#library(pROC)
library(cowplot)
library(survival)
ukb <- fread('E:\\English Papers for writing\\have done\\ALCO-MAFLD\\ukb_alco_mafld.csv')
ukb2 <- ukb %>% mutate_if(
  is.character,
  function(x) ifelse(x == '',NA,x)) %>% 
  mutate(
    # DiaInterval = time_length(ymd(f.41280.0.0)-ymd(f.53.0.0),'day'),
    DeathInterval = time_length(ymd(f.40000.0.0)-ymd(f.53.0.0),'day'),
    CancerInterval = time_length(ymd(f.40005.0.0)-ymd(f.53.0.0),'day'),
    survive = ifelse(is.na(f.40000.0.0),'survive','die')
  )
ukb3 <- ukb2 %>% rename(
  eid = f.eid,sex = f.31.0.0,YearOfBirth = f.34.0.0,
  attendingCentre = f.54.0.0,
  MonthOfBirth = f.52.0.0,
  WC = f.30000.0.0,RC = f.30010.0.0,HB = f.30020.0.0,
  HaePer = f.30030.0.0,MCV = f.30040.0.0,MCH = f.30050.0.0,
  MCHC = f.30060.0.0,RBCDW = f.30070.0.0,PLT = f.30080.0.0,
  PLTC = f.30090.0.0,MPV = f.30100.0.0,PLTDW = f.30100.0.0,
  LYM = f.30120.0.0,MON = f.30130.0.0,NEU = f.30140.0.0,
  EOS = f.30150.0.0,BAS = f.30160.0.0,NRBC = f.30170.0.0,
  LYMPer = f.30180.0.0,MonPer = f.30190.0.0,NeuPer = f.30200.0.0,
  EsoPer = f.30210.0.0,BasPer = f.30220.0.0,NRBCPer = f.30230.0.0,
  RetPer = f.30240.0.0,RET = f.30250.0.0,MRV = f.30260.0.0,
  MSCV = f.30270.0.0,IRF = f.30280.0.0,HLSRPer = f.30290.0.0,
  HLSR = f.30300.0.0,MALBUR = f.30500.0.0,CREUR = f.30510.0.0,
  KaUR = f.30520.0.0,NaUR = f.30530.0.0,ALB = f.30600.0.0,
  ALP = f.30610.0.0,ALT = f.30620.0.0,ApoA = f.30630.0.0,
  ApoB = f.30640.0.0,AST = f.30650.0.0,DBi = f.30660.0.0,
  UREA = f.30670.0.0,Ca = f.30680.0.0,CHL = f.30690.0.0,
  CRE = f.30700.0.0,CRP = f.30710.0.0,CysC = f.30720.0.0,
  GGT = f.30730.0.0,Glucose = f.30740.0.0,HbA1c = f.30750.0.0,
  HDL = f.30760.0.0,IGF1 = f.30770.0.0,LDL = f.30780.0.0,
  LipA = f.30790.0.0,OES = f.30800.0.0,PHOS = f.30810.0.0,
  RHE = f.30820.0.0,SHBG = f.30830.0.0,TBi = f.30840.0.0,
  TES = f.30850.0.0,TPro = f.30860.0.0,TRG = f.30870.0.0,
  Urate = f.30880.0.0,VD = f.30890.0.0,
  DateDeath = f.40000.0.0,DeathCause = f.40001.0.0,
  DateCaDia = f.40005.0.0,TypeofCa = f.40006.0.0,
  AgeDeath = f.40007.0.0,AgeCaDia = f.40008.0.0,
  HistoCa = f.40011.0.0
)
ukb3 <- ukb3 %>% rename(handGripStrength = f.47.0.0,waistCir = f.48.0.0,hipCir = f.49.0.0,
                        standingHeight = f.50.0.0,seatedHeight = f.51.0.0,dateAttendingCentre = f.53.0.0,
                        cancerAgeOccurred = f.84.0.0,noncancerAgeOccurred = f.87.0.0,
                        BirthWeightType = f.120.0.0,no.selfreportCA = f.134.0.0,
                        no.selfreportNonCA = f.135.0.0,no.selfreportTreat = f.137.0.0,
                        dateLossFollowup = f.191.0.0, averageTotalIncome = f.738.0.0,
                        ageCompleteFulltimeEdu = f.845.0.0,
                        no.daysWalk10 = f.864.0.0,
                        no.daysModeAct = f.884.0.0,
                        no.daysVigorAct = f.904.0.0,
                        sleepDuration = f.1160.0.0,insomnia = f.1200.0.0,
                        currentSmoke = f.1239.0.0,pastSmoke = f.1249.0.0,
                        aolIntakeFreq = f.1558.0.0,redwine1 = f.1568.0.0,
                        whitewine1 = f.1578.0.0,beer1 = f.1588.0.0,
                        spirits1 = f.1598.0.0,fortifiedwine1 = f.1608.0.0,
                        otherwine1 = f.5364.0.0,
                        redwine2 = f.4407.0.0,
                        whitewine2 = f.4418.0.0,beer2 = f.4429.0.0,
                        spirits2 = f.4440.0.0,fortifiedwine2 = f.4451.0.0,
                        otherwine2 = f.4462.0.0,
                        wineChange = f.1628.0.0,
                        CompBodysizeAt10 = f.1687.0.0,CompHeightsizeAt10= f.1697.0.0,
                        overallHealth = f.2178.0.0,longstandingIll = f.2188.0.0,
                        weightChange = f.2306.0.0,hadMajorOperation = f.2415.0.0,
                        diabetesBydoctor = f.2443.0.0,caBydoctor = f.2453.0.0,
                        otherSeriousIllBydoc = f.2473.0.0,
                        hadPrescription = f.2492.0.0,ageMenarche = f.2714.0.0,
                        hadOralContraceptivePill = f.2784.0.0,ageContraceptivePill=f.2794.0.0,
                        ageLastUseCP = f.2804.0.0,ageSmoking = f.2867.0.0,
                        forcedVitalCapacity = f.3062.0.0,FEV1 = f.3063.0.0,
                        PeakExpiratoryFlow = f.3064.0.0,coffeeDrink = f.3089.0.0,
                        DBP = f.4079.0.0,SBP = f.4080.0.0,
                        Nevereateggsdairywheatsugar = f.6144.0.0,
                        allowance = f.6146.0.0,dentalProbelm = f.6149.0.0,
                        vascularDisBydoc = f.6150.0.0,otherDisBydoc = f.6152.0.0,
                        medCholBPDiaF = f.6153.0.0,medPain = f.6154.0.0,
                        vitaMineralSupp = f.6155.0.0,medCholBPDiaM = f.6177.0.0,
                        otherDietarySupp = f.6179.0.0,birthweight = f.20022.0.0,
                        illnessMom = f.20110.0.0,smokingStatus = f.20116.0.0,
                        alcoholStatus = f.20117.0.0, NeuroticismScore = f.20127.0.0,
                        ethnic = f.21000.0.0,BMI = f.21001.0.0,weight = f.21002.0.0,
                        ageAssessment = f.21003.0.0,ageRecruitment = f.21022.0.0,
                        IBS = f.21024.0.0,
                        bodyFatPerc = f.23099.0.0,wholeBodyFatMass = f.23100.0.0,
                        wholeBodyFatFreeMass = f.23101.0.0, BMI2 = f.23104.0.0,
                        multipleDeprivationIndex = f.26410.0.0,
                        IncomeScore = f.26411.0.0,
                        healthScore = f.26413.0.0,
                        educationScore = f.26414.0.0,
                        housingScore = f.26415.0.0,
                        crimeScore = f.26416.0.0,
                        MPV = PLTDW, PLTDW = f.30110.0.0)

grs <- fread('snps_query_data.csv')
pdff <- fread('data_for_PDFF_GWAS.csv')

###explained proportions by snps
pdff <- pdff %>% left_join(grs,by = c('f.eid'='ID_1'))
summary(lm(FRPDFF~rs738408+rs58542926+rs429358+rs188247550+rs28601761+rs7096937+rs1260326+rs2642438+rs1229984+rs79905393+rs62226381,
           data = pdff))
##R2 = 4.8%
pdff <- pdff %>% mutate(
  HFCGRS=rs58542926*0.3080+rs188247550*0.2825+rs738408*0.2375+rs1229984*0.1673+rs2642438*0.0604+rs62226381*0.0547
  -rs1260326*0.0586-rs28601761*0.0618-rs7096937*0.0654-rs79905393*0.1130-rs429358*0.1150
)
cor.test(pdff$FRPDFF,pdff$HFCGRS)
ggplot(pdff,aes(HFCGRS,FRPDFF))+
  geom_jitter()+
  geom_smooth(method = "lm")+
  annotate('text',x = 1,y=-3.5,label='rho = 0.219, p < 2.2e-16',size=6)

grs <- grs %>% filter(!is.na(rs738408))
ukb4 <- ukb3 %>% left_join(grs,by = c('eid'='ID_1'))
ukb5 <-  ukb4 %>% filter(!eid %in% pdff$f.eid)
ukb6 <- ukb5 %>% filter(!is.na(rs738408))
ukb6 <- ukb6 %>%
  mutate(
    age = year(dateAttendingCentre)-YearOfBirth,
    education = ifelse(!is.na(educationScore),educationScore,
                       ifelse(!is.na(f.26421.0.0),f.26421.0.0,f.26431.0.0))) %>% 
  mutate(
    education = ifelse(is.na(education),median(education,na.rm = T),education),
    averageTotalIncome = ifelse(is.na(averageTotalIncome)|averageTotalIncome<0,
                                median(averageTotalIncome,na.rm = T),averageTotalIncome),
    aolIntakeFreq = ifelse(!aolIntakeFreq %in% 1:6,3,aolIntakeFreq),
    BMI = ifelse(is.na(BMI),median(BMI,na.rm = T),BMI),
    smokingStatus = ifelse(!smokingStatus%in% 0:2,0,smokingStatus),
    no.daysModeAct = ifelse(no.daysModeAct<0|is.na(no.daysModeAct),3,no.daysModeAct),
    HbA1c = ifelse(is.na(HbA1c),median(HbA1c,na.rm = T),HbA1c),
    Glucose = ifelse(is.na(Glucose),median(Glucose,na.rm = T),Glucose),
    CHL = ifelse(is.na(CHL),median(CHL,na.rm = T),CHL),
    LDL = ifelse(is.na(LDL),median(LDL,na.rm = T),LDL),
    HDL = ifelse(is.na(HDL),median(HDL,na.rm = T),HDL),
    ALT = ifelse(is.na(ALT),median(ALT,na.rm = T),ALT),
    AST = ifelse(is.na(AST),median(AST,na.rm = T),AST),
    SBP = ifelse(is.na(SBP),median(SBP,na.rm = T),SBP),
    DBP = ifelse(is.na(DBP),median(DBP,na.rm = T),DBP),
    GGT = ifelse(is.na(GGT),median(GGT,na.rm = T),GGT),
    TRG = ifelse(is.na(TRG),median(TRG,na.rm = T),TRG),
    PLT = ifelse(is.na(PLT),median(PLT,na.rm = T),PLT),
    ALB = ifelse(is.na(ALB),median(ALB,na.rm = T),ALB),
    waistCir = ifelse(is.na(waistCir),median(waistCir,na.rm = T),waistCir)
  ) %>% mutate(
    no.daysModeAct = ifelse(no.daysModeAct<2,'0-1',ifelse(no.daysModeAct<5,'2-4','5-7')),
    aolIntakeFreq = ifelse(aolIntakeFreq>=5,'Never',ifelse(aolIntakeFreq<=2,'Excess','Moderate')),
    aolIntakeFreq = factor(aolIntakeFreq,levels = c('Never','Moderate','Excess')),
    TRG2 = TRG/0.0113,
    FLI = exp(0.953*log(TRG2)+0.139*BMI+0.718*log(GGT)+
                0.053*waistCir-15.745)/(1+exp(0.953*log(TRG2)+0.139*BMI+0.718*log(GGT)+
                                                0.053*waistCir-15.745))*100,
    HFCGRS=rs58542926*0.3080+rs188247550*0.2825+rs738408*0.2375+rs1229984*0.1673+rs2642438*0.0604+rs62226381*0.0547
    -rs1260326*0.0586-rs28601761*0.0618-rs7096937*0.0654-rs79905393*0.1130-rs429358*0.1150
  )
ukb7 <- ukb6 %>% filter(ethnic %in% c(1001,1002,1003))
ggplot(ukb7,aes(HFCGRS))+
  geom_density(adjust=3,color='#00468b',linewidth=2)+
  theme_bw()+
  labs(x = 'PRS of hepatic fat content')
ukb7 <- ukb7 %>% rename(sex = sex.x)
ukb7 <- ukb7 %>% mutate(
  loweringCHL = ifelse(((sex==1 & medCholBPDiaM==1) |(sex==0 & medCholBPDiaF==1)),1,0),
  loweringINS = ifelse(((sex==1 & medCholBPDiaM==3) |(sex==0 & medCholBPDiaF==3)),1,0),
  loweringBP  = ifelse(((sex==1 & medCholBPDiaM==2) |(sex==0 & medCholBPDiaF==2)),1,0),
  VD_supple =   ifelse(vitaMineralSupp==4,1,0)
)
fwrite(ukb7,'analyzingData.csv')  
##disease phenotypes
ukb7 <- fread('analyzingData.csv')
phenos <- fread('phenos.csv')
phenos2 <- phenos %>% filter(eid %in% ukb7$eid)
phenoGroup <- read.csv('phewas_results_old.csv')
phenos2 <- phenos2 %>% dplyr::left_join(phenoGroup %>% dplyr::select(PheCode,PhenotypeGroup),by='PheCode')

phenos_freq <- as.data.frame(table(phenos2$PheCode)) %>% 
  filter(Freq>=20)

phenos_freq <- phenos_freq %>% mutate(Var1=as.numeric(as.character(Var1))) %>%
  rename(PheCode = Var1)
phenos3 <- phenos2 %>% filter(PheCode %in% phenos_freq$PheCode)
PheCodes <- unique(phenos3$PheCode)
Phenotypes <- unique(phenos3$Phenotype)
tempData <- ukb7 %>% dplyr::select(eid,HFCGRS,sex,age,smokingStatus,averageTotalIncome,aolIntakeFreq,
                            no.daysModeAct,education,BMI)
tempData2 <- tempData %>% dplyr::left_join(phenos3 %>% dplyr::select(eid,PheCode,PhenotypeGroup),by='eid')

df2 <- data.frame()
system.time(
  for(i in 1:length(PheCodes)){
    case <- tempData2 %>% filter(PheCode == PheCodes[i]) %>% mutate(pheno =1)
    phenotypeGroup <- tempData2 %>% filter(PheCode== PheCodes[i]) %>% select(PhenotypeGroup)
    control <- tempData2 %>% filter(PhenotypeGroup != phenotypeGroup$PhenotypeGroup[1]) %>% 
      filter(!eid %in% case$eid) %>%
      distinct(eid,.keep_all = TRUE) %>% 
      mutate(pheno=0)
    dat1 <- case %>% bind_rows(control)
    fit <- glm(pheno~scale(HFCGRS)+sex+age+smokingStatus+averageTotalIncome+aolIntakeFreq+
                 no.daysModeAct+education+BMI,
               family = binomial,data = dat1)
    fit_estimate <- broom::tidy(fit)
    df <- data.frame(PheCode = PheCodes[i],
                     phenotype=Phenotypes[i],
                     n_case = dim(case)[1],
                     n_control = dim(control)[1],
                     beta = fit_estimate$estimate[2],
                     se=fit_estimate$std.error[2],
                     p = fit_estimate$p.value[2])
    df2 <- bind_rows(df2,df)
    print(paste(i,' was finished!!!'))
  }
)

phewas_res <-df2 %>% left_join(phenoGroup%>%select(PheCode,PhenotypeGroup),
                               by = 'PheCode') 
phewas_res <- phewas_res %>% 
  mutate(
    q = p.adjust(p,'fdr'),  
    logp = -log10(q),
    PhenotypeGroup = factor(PhenotypeGroup,
                            levels = unique(phenoGroup$PhenotypeGroup))
  )

fwrite(phewas_res,'phewas_results.csv',row.names = F)

phewas_res <- fread('phewas_results.csv')
phewas_res <- phewas_res %>% 
  mutate(
    PhenotypeGroup = factor(PhenotypeGroup,
                            levels = unique(phenoGroup$PhenotypeGroup))
  )
#####phewas plot
library(RColorBrewer)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, 
                           rownames(qual_col_pals))) 

ggplot()+
  geom_hline(yintercept = 1.3,color = 'red',linetype=2)+
  geom_jitter(data=phewas_res%>% filter(q>=0.05),
              aes(PhenotypeGroup,logp,color = PhenotypeGroup),
              size=2)+
  geom_jitter(data=phewas_res%>% filter(q<0.05,logp<15),
              aes(PhenotypeGroup,logp,color = PhenotypeGroup),
              size=3)+
  theme_classic()+
  ggrepel::geom_text_repel(data = phewas_res %>% filter(q<0.05,logp<15),
                           aes(PhenotypeGroup,logp,label=phenotype),max.overlaps = 20,
                           angle=30)+
  scale_x_discrete(expand = expand_scale(mult=0.01,add = c(0.5,0.5)))+
  scale_y_continuous(breaks = seq(0,14,by=2),
                     expand = expand_scale(mult = 0.02),
                     limits = c(0,15))+
  labs(x = NULL,
       y = '-log10(FDR)')+
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_text(color='black',angle = 45,size=12,hjust=1),
        axis.text.y = element_text(color='black',size=12),
        axis.title = element_text(size=13),
        legend.position = 'none')+
  scale_color_manual(values = c(col_vector[10:23]))+
  guides(color = guide_legend(title=NULL,ncol = 4))

ggplot()+
  geom_jitter(data=phewas_res%>% filter(q>=0.05),
              aes(PhenotypeGroup,logp,color = PhenotypeGroup),
              size=2)+
  geom_jitter(data=phewas_res%>% filter(q<0.05,logp>=15),
              aes(PhenotypeGroup,logp,color = PhenotypeGroup),
              size=3)+
  theme_classic()+
  ggrepel::geom_text_repel(data = phewas_res %>% filter(q<0.05,logp>=15),
                           aes(PhenotypeGroup,logp,label=phenotype),max.overlaps = 20,
                           angle=30)+
  scale_x_discrete(expand = expand_scale(mult=0.01,add = c(0.5,0.5)))+
  scale_y_continuous(breaks = seq(15,65,by=10),
                     expand = expand_scale(mult = 0.02),
                     limits = c(15,65))+
  labs(x = NULL,
       y = '-log10(FDR)')+
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(color='black',size=12),
        axis.title = element_text(size=13),
        legend.position = 'none')+
  scale_color_manual(values = c(col_vector[10:23]))+
  guides(color = guide_legend(title=NULL,ncol = 4))


####phenotypes OR plot
phewas_res <- read.csv('phewas_results.csv')
or <- phewas_res %>% filter(q<0.05) %>%
  mutate(
    OR = exp(beta),
    lwr= exp(beta-1.96*se),
    upr= exp(beta+1.96*se),
    phenotype= factor(phenotype,levels = phenotype[order(OR,decreasing = T)])
  )
or <- or %>% mutate(
  `OR (95% CI)` = sprintf("%.2f (%.2f to %.2f)",
                          or$OR, or$lwr, or$upr)
)

a <- ggplot(or,aes(OR,phenotype))+
  geom_vline(xintercept = 1,linetype=2)+
  geom_errorbarh(aes(xmin=lwr,xmax=upr),height=0)+
  geom_point(size=2.5)+
  scale_x_continuous(breaks = seq(0.7,1.8,by=.1),
                     limits = c(0.7,2.5))+
  scale_y_discrete(expand = expansion(add = c(0.5,1.5)))+
  #scale_color_manual(name = 'Phenotype group',values = col_vector[c(10,11,12,14,16,18,19,21,22)])+
  theme_test()+
  labs(x = 'Odds ratio (95% CI)',y=NULL)+
  geom_text(aes(x = 2.,y=phenotype,label = `OR (95% CI)`),color='black')+
  geom_text(aes(x=2.45,y=phenotype,label=scales::scientific(q),fontface = 'italic'),color='black')+
  annotate('text',x = 2,y=29,label = 'Odds ratio (95% CI)',color='black',fontface='bold')+
  annotate('text',x = 2.45,y=29,label = 'FDR',color='black',fontface='bold')+
  theme(axis.text = element_text(color='black',size = 11),
        legend.text = element_text(size=11),
        axis.title = element_text(size=13))


###odds ratio for biomarkers
betas <- read.csv('biomarkers_HFCGRS.csv')
betas2 <- betas %>% filter(fdr<0.05) %>%
  mutate(
    BiomarkerNames= factor(BiomarkerNames,levels = BiomarkerNames[order(beta,decreasing = T)]),
    lwr = beta-1.96*se,
    upr = beta+1.96*se
  )
betas2 <- betas2 %>% mutate(
  `Beta (95% CI)` = sprintf("%.3f (%.3f to %.3f)",
                          betas2$beta, betas2$lwr, betas2$upr)
)

b <- ggplot(betas2,aes(beta,BiomarkerNames))+
  geom_vline(xintercept = 0,linetype=2)+
  geom_errorbarh(aes(xmin=lwr,xmax=upr),height=0)+
  geom_point(size=2.5)+
  scale_x_continuous(breaks = seq(-0.04,0.09,by=.02),
                     limits = c(-0.04,0.18))+
  scale_y_discrete(expand = expansion(add = c(0.5,1.5)))+
  theme_test()+
  labs(x = 'Beta (95% CI)',y=NULL)+
  geom_text(aes(x = 0.12,y=BiomarkerNames,label = `Beta (95% CI)`),color='black')+
  geom_text(aes(x=0.175,y=BiomarkerNames,label=scales::scientific(fdr),fontface = 'italic'),color='black')+
  annotate('text',x = 0.12,y=43,label = 'Beta (95% CI)',color='black',fontface='bold')+
  annotate('text',x = 0.175,y=43,label = 'FDR',color='black',fontface='bold')+
  theme(axis.text = element_text(color='black',size = 11),
        legend.text = element_text(size=11),
        axis.title = element_text(size=13))
plot_grid(a,b,ncol = 1,labels = c('B','C'),align = "v",rel_heights = c(.8,1))
####TWO-SAMPLE MR--------------------------------------
library(TwoSampleMR)
library(MRPRESSO)

mr_scatter_plot2 <- function (mr_results, dat) 
{
  requireNamespace("ggplot2", quietly = TRUE)
  requireNamespace("plyr", quietly = TRUE)
  mrres <- plyr::dlply(dat, c("id.exposure", "id.outcome"), 
                       function(d) {
                         d <- plyr::mutate(d)
                         if (nrow(d) < 2 | sum(d$mr_keep) == 0) {
                           return(blank_plot("Insufficient number of SNPs"))
                         }
                         d <- subset(d, mr_keep)
                         index <- d$beta.exposure < 0
                         d$beta.exposure[index] <- d$beta.exposure[index] * 
                           -1
                         d$beta.outcome[index] <- d$beta.outcome[index] * 
                           -1
                         mrres <- subset(mr_results, id.exposure == d$id.exposure[1] & 
                                           id.outcome == d$id.outcome[1])
                         mrres$a <- 0
                         if ("MR Egger" %in% mrres$method) {
                           temp <- mr_egger_regression(d$beta.exposure, 
                                                       d$beta.outcome, d$se.exposure, d$se.outcome, 
                                                       default_parameters())
                           mrres$a[mrres$method == "MR Egger"] <- temp$b_i
                         }
                         if ("MR Egger (bootstrap)" %in% mrres$method) {
                           temp <- mr_egger_regression_bootstrap(d$beta.exposure, 
                                                                 d$beta.outcome, d$se.exposure, d$se.outcome, 
                                                                 default_parameters())
                           mrres$a[mrres$method == "MR Egger (bootstrap)"] <- temp$b_i
                         }
                         ggplot2::ggplot(data = d, ggplot2::aes(x = beta.exposure, 
                                                                y = beta.outcome)) + 
                           ggplot2::geom_errorbar(ggplot2::aes(ymin = beta.outcome - se.outcome, 
                                                               ymax = beta.outcome + se.outcome), 
                                                  colour = "grey", width = 0) + 
                           ggplot2::geom_errorbarh(ggplot2::aes(xmin = beta.exposure -se.exposure, 
                                                                xmax = beta.exposure + se.exposure), 
                                                   colour = "grey", height = 0) + 
                           ggplot2::geom_point() + 
                           ggplot2::geom_abline(data = mrres, ggplot2::aes(intercept = a, 
                                                                           slope = b, colour = method), 
                                                show.legend = TRUE,size=1.2) + 
                           ggsci::scale_color_npg() + 
                           ggplot2::labs(colour = "MR Test", 
                                         x = paste("SNP effect on exposure"), 
                                         y = paste("SNP effect on outcome")) + 
                           ggplot2::theme_test()+
                           ggplot2::theme(legend.position = "top", 
                                          legend.direction = "vertical") + 
                           ggplot2::guides(colour = ggplot2::guide_legend(ncol = 2))
                       })
  mrres
}
mr_scatter_plot3 <- function (mr_results, dat) 
{
  requireNamespace("ggplot2", quietly = TRUE)
  requireNamespace("plyr", quietly = TRUE)
  mrres <- plyr::dlply(dat, c("id.exposure", "id.outcome"), 
                       function(d) {
                         d <- plyr::mutate(d)
                         if (nrow(d) < 2 | sum(d$mr_keep) == 0) {
                           return(blank_plot("Insufficient number of SNPs"))
                         }
                         d <- subset(d, mr_keep)
                         index <- d$beta.exposure < 0
                         d$beta.exposure[index] <- d$beta.exposure[index] * 
                           -1
                         d$beta.outcome[index] <- d$beta.outcome[index] * 
                           -1
                         mrres <- subset(mr_results, id.exposure == d$id.exposure[1] & 
                                           id.outcome == d$id.outcome[1])
                         mrres$a <- 0
                         if ("MR Egger" %in% mrres$method) {
                           temp <- mr_egger_regression(d$beta.exposure, 
                                                       d$beta.outcome, d$se.exposure, d$se.outcome, 
                                                       default_parameters())
                           mrres$a[mrres$method == "MR Egger"] <- temp$b_i
                         }
                         if ("MR Egger (bootstrap)" %in% mrres$method) {
                           temp <- mr_egger_regression_bootstrap(d$beta.exposure, 
                                                                 d$beta.outcome, d$se.exposure, d$se.outcome, 
                                                                 default_parameters())
                           mrres$a[mrres$method == "MR Egger (bootstrap)"] <- temp$b_i
                         }
                         ggplot2::ggplot(data = d, ggplot2::aes(x = beta.exposure, 
                                                                y = beta.outcome)) + 
                           ggplot2::geom_errorbar(ggplot2::aes(ymin = beta.outcome - se.outcome, 
                                                               ymax = beta.outcome + se.outcome), 
                                                  colour = "grey", width = 0) + 
                           ggplot2::geom_errorbarh(ggplot2::aes(xmin = beta.exposure -se.exposure, 
                                                                xmax = beta.exposure + se.exposure), 
                                                   colour = "grey", height = 0) + 
                           ggplot2::geom_point() + 
                           ggplot2::geom_abline(data = mrres, ggplot2::aes(intercept = a, 
                                                                           slope = b, colour = method), 
                                                show.legend = TRUE,size=1.2) + 
                           ggsci::scale_color_npg() + 
                           ggplot2::labs(colour = "MR Test", 
                                         x = paste("SNP effect on exposure"), 
                                         y = paste("SNP effect on outcome")) + 
                           ggplot2::theme_test()+
                           ggplot2::theme(legend.position = "none") + 
                           ggplot2::guides(colour = ggplot2::guide_legend(ncol = 2))
                       })
  mrres
}


hfc <- read.csv('independent_sig_snps.csv')
hfc_mr <- format_data(hfc,
                      type = 'exposure',
                      snp_col =  'snp',
                      beta_col = 'beta',
                      se_col = 'se',
                      eaf_col = 'EAF',
                      effect_allele_col = 'Alt_allele',
                      other_allele_col = 'Ref_allele',
                      pval_col = 'pval',
                      chr_col = 'chr',
                      pos_col = 'pos')

r2 = 2*hfc_mr$beta.exposure^2*hfc_mr$eaf.exposure*(1-hfc_mr$eaf.exposure)
mean(r2*(31377-2)/(1-r2))
#=136.9
###hepatitis nos
CHP <- fread('finngen_R7_K11_CHRONHEP.gz')
CHP2 <- format_data(CHP,
                    type = 'outcome',
                    snps = hfc_mr$snp,
                    snp_col =  'rsids',
                    beta_col = 'beta',
                    se_col = 'sebeta',
                    eaf_col = 'af_alt',
                    effect_allele_col = 'alt',
                    other_allele_col = 'ref',
                    pval_col = 'pval',
                    chr_col = '#chrom',
                    pos_col = 'pos')
dat <- harmonise_data(hfc_mr,CHP2)
write.csv(dat,'chp.csv',row.names = F)
dat <- read.csv('chp.csv')
dat <- dat %>% mutate(mr_keep=T)
mr_res = mr(dat,method_list = c('mr_ivw_fe',
                                'mr_egger_regression',
                                'mr_weighted_median',
                                'mr_weighted_mode'))
mr_heterogeneity(dat)
#6.5 0.167
mr_pleiotropy_test(dat)
df <- data.frame(
  exposure  = 'MRI-HFC',
  outcome   = 'Hepatitis NOS',
  method    = mr_res$method,
  nsnp      = mr_res$nsnp,
  b         = mr_res$b,
  se        = mr_res$se,
  pval      = mr_res$pval
)

a <- mr_scatter_plot3(mr_res,dat)
a2 <- a$LKFN1o.q04Z5v+labs(title = 'Hepatitis NOS')

rm(list = c('CHP','CHP2'))

##liver cancer
LC <- fread('finngen_R7_C3_LIVER_INTRAHEPATIC_BILE_DUCTS_EXALLC.gz')
LC2 <- format_data(LC,
                   type = 'outcome',
                   snps = hfc_mr$snp,
                   snp_col =  'rsids',
                   beta_col = 'beta',
                   se_col = 'sebeta',
                   eaf_col = 'af_alt',
                   effect_allele_col = 'alt',
                   other_allele_col = 'ref',
                   pval_col = 'pval',
                   chr_col = '#chrom',
                   pos_col = 'pos')
dat <- harmonise_data(hfc_mr,LC2)
write.csv(dat,'mr_liverCancer.csv',row.names = F)
dat <- dat %>% mutate(mr_keep=T)
mr_res = mr(dat,method_list = c('mr_ivw_fe',
                                'mr_egger_regression',
                                'mr_weighted_median',
                                'mr_weighted_mode'))
mr_heterogeneity(dat)
#6.5 0.167
mr_pleiotropy_test(dat)
df <- df %>% bind_rows(
  data.frame(
  exposure  = 'MRI-HFC',
  outcome   = 'Primary liver cancer',
  method    = mr_res$method,
  nsnp      = mr_res$nsnp,
  b         = mr_res$b,
  se        = mr_res$se,
  pval      = mr_res$pval
))

b <- mr_scatter_plot3(mr_res,dat)
b2 <- b$LKFN1o.MMxQIo+labs(title = 'Primary liver cancer')

rm(list = c('LC','LC2'))

##Polycythemia vera
PV <- fread('finngen_R7_POLYCYTVERA.gz')
PV2 <- format_data(PV,
                   type = 'outcome',
                   snps = hfc_mr$snp,
                   snp_col =  'rsids',
                   beta_col = 'beta',
                   se_col = 'sebeta',
                   eaf_col = 'af_alt',
                   effect_allele_col = 'alt',
                   other_allele_col = 'ref',
                   pval_col = 'pval',
                   chr_col = '#chrom',
                   pos_col = 'pos')
dat <- harmonise_data(hfc_mr,PV2)
write.csv(dat,'mr_Polycythemia vera.csv',row.names = F)
dat <- dat %>% mutate(mr_keep=T)
mr_res = mr(dat,method_list = c('mr_ivw_fe',
                                'mr_egger_regression',
                                'mr_weighted_median',
                                'mr_weighted_mode'))
mr_heterogeneity(dat)
mr_pleiotropy_test(dat)
df <- df %>% bind_rows(
  data.frame(
    exposure  = 'MRI-HFC',
    outcome   = 'Polycythemia vera',
    method    = mr_res$method,
    nsnp      = mr_res$nsnp,
    b         = mr_res$b,
    se        = mr_res$se,
    pval      = mr_res$pval
  ))

# b <- mr_scatter_plot3(mr_res,dat)
# b2 <- b$LKFN1o.MMxQIo+labs(title = 'Polycythemia vera')

rm(list = c('PV','PV2'))

###T2D
t2d <- fread('Mahajan.NatGenet2018b.T2D-noUKBB.European.txt')
hfc <- hfc %>% mutate(SNP = paste0(chr,":",pos))
t2d2 <- t2d %>% filter(SNP %in% hfc$SNP)
t2d2 <- t2d2 %>% mutate(snp = hfc$snp)
t2d2 <- format_data(t2d2,
                    type = 'outcome',
                    snp_col =  'snp',
                    beta_col = 'Beta',
                    se_col = 'SE',
                    eaf_col = 'EAF',
                    effect_allele_col = 'EA',
                    other_allele_col = 'NEA',
                    pval_col = 'Pvalue',
                    chr_col = 'Chr',
                    pos_col = 'Pos')
dat <- harmonise_data(hfc_mr,t2d2)
write.csv(dat,'mr_T2D.csv',row.names = F)
dat <- dat %>% mutate(mr_keep=T)
mr_res = mr(dat,method_list = c('mr_ivw_fe',
                                'mr_egger_regression',
                                'mr_weighted_median',
                                'mr_weighted_mode'))
mr_heterogeneity(dat)
mr_pleiotropy_test(dat)
df <- df %>% bind_rows(
  data.frame(
    exposure  = 'MRI-HFC',
    outcome   = 'Type 2 diabetes',
    method    = mr_res$method,
    nsnp      = mr_res$nsnp,
    b         = mr_res$b,
    se        = mr_res$se,
    pval      = mr_res$pval
  ))
cc <- mr_scatter_plot3(mr_res,dat)
cc2 <- cc$LKFN1o.NjrTVq+labs(title = 'Type 2 diabetes')
rm(list = c('t2d','t2d2'))

###Hypercholesterolemia
HLA <- fread('DYSLIPID_GERA.txt.gz')
HLA2 <- format_data(HLA,
                    type = 'outcome',
                    snps = hfc_mr$snp,
                    snp_col =  'name',
                    beta_col = 'frequentist_add_beta_1',
                    se_col = 'frequentist_add_se_1',
                    eaf_col = 'all_maf',
                    effect_allele_col = 'alleleB',
                    other_allele_col = 'alleleA',
                    pval_col = 'frequentist_add_pvalue',
                    chr_col = 'chr',
                    pos_col = 'position')
dat <- harmonise_data(hfc_mr,HLA2)
write.csv(dat,'mr_Hypercholesterolemia.csv',row.names = F)
dat <- dat %>% mutate(mr_keep=T)
mr_res = mr(dat,method_list = c('mr_ivw_fe',
                                'mr_egger_regression',
                                'mr_weighted_median',
                                'mr_weighted_mode'))
mr_heterogeneity(dat)
mr_pleiotropy_test(dat)
df <- df %>% bind_rows(
  data.frame(
    exposure  = 'MRI-HFC',
    outcome   = 'Hypercholesterolemia',
    method    = mr_res$method,
    nsnp      = mr_res$nsnp,
    b         = mr_res$b,
    se        = mr_res$se,
    pval      = mr_res$pval
  ))
d <- mr_scatter_plot3(mr_res,dat)
d2 <- d$LKFN1o.4TmsH2+labs(title = 'Hypercholesterolemia')
rm(list = c('HLA','HLA2'))

###dementia
DEM <- fread('finngen_R7_F5_DEMENTIA.gz')
DEM2 <- format_data(DEM,
                    type = 'outcome',
                    snps = hfc_mr$snp,
                    snp_col =  'rsids',
                    beta_col = 'beta',
                    se_col = 'sebeta',
                    eaf_col = 'af_alt',
                    effect_allele_col = 'alt',
                    other_allele_col = 'ref',
                    pval_col = 'pval',
                    chr_col = '#chrom',
                    pos_col = 'pos')
dat <- harmonise_data(hfc_mr,DEM2)
write.csv(dat,'mr_dementia.csv',row.names = F)
dat <- dat %>% mutate(mr_keep=T)
mr_res = mr(dat,method_list = c('mr_ivw_fe',
                                'mr_egger_regression',
                                'mr_weighted_median',
                                'mr_weighted_mode'))
mr_heterogeneity(dat)
mr_pleiotropy_test(dat)
df <- df %>% bind_rows(
  data.frame(
    exposure  = 'MRI-HFC',
    outcome   = 'Dementia',
    method    = mr_res$method,
    nsnp      = mr_res$nsnp,
    b         = mr_res$b,
    se        = mr_res$se,
    pval      = mr_res$pval
  ))
# d <- mr_scatter_plot3(mr_res,dat)
# d2 <- d$LKFN1o.4TmsH2+labs(title = 'Hypercholesterolemia')
rm(list = c('DEM','DEM2'))

###AD
AD <- fread('PGCALZ2ExcludingUKBand23andME_METALInverseVariance_MetaAnalysis.txt.gz')
AD2 <- AD %>% mutate(snp=paste0(chromosome,":",base_pair_location)) %>%
  filter(snp %in% hfc$SNP)
AD2 <- AD2 %>% mutate(SNP = hfc$snp)
AD2 <- format_data(AD2,
                    type = 'outcome',
                    snp_col =  'SNP',
                    beta_col = 'beta',
                    se_col = 'standard_error',
                    eaf_col = 'effect_allele_frequency',
                    effect_allele_col = 'effect_allele',
                    other_allele_col = 'other_allele',
                    pval_col = 'p_value',
                    chr_col = 'chromosome',
                    pos_col = 'base_pair_location')
dat <- harmonise_data(hfc_mr,AD2)
write.csv(dat,'mr_AD.csv',row.names = F)
dat <- dat %>% mutate(mr_keep=T)
mr_res = mr(dat,method_list = c('mr_ivw_fe',
                                'mr_egger_regression',
                                'mr_weighted_median',
                                'mr_weighted_mode'))
mr_heterogeneity(dat)
mr_pleiotropy_test(dat)
df <- df %>% bind_rows(
  data.frame(
    exposure  = 'MRI-HFC',
    outcome   = "Alzheimer’s disease",
    method    = mr_res$method,
    nsnp      = mr_res$nsnp,
    b         = mr_res$b,
    se        = mr_res$se,
    pval      = mr_res$pval
  ))
# d <- mr_scatter_plot3(mr_res,dat)
# d2 <- d$LKFN1o.4TmsH2+labs(title = 'Hypercholesterolemia')
rm(list = c('AD','AD2'))

###vascular dementia
VDE <- fread('finngen_R7_F5_VASCDEM.gz')
VDE2 <- format_data(VDE,
                    type = 'outcome',
                    snps = hfc_mr$snp,
                    snp_col =  'rsids',
                    beta_col = 'beta',
                    se_col = 'sebeta',
                    eaf_col = 'af_alt',
                    effect_allele_col = 'alt',
                    other_allele_col = 'ref',
                    pval_col = 'pval',
                    chr_col = '#chrom',
                    pos_col = 'pos')
dat <- harmonise_data(hfc_mr,VDE2)
write.csv(dat,'mr_vascular dementia.csv',row.names = F)
dat <- dat %>% mutate(mr_keep=T)
mr_res = mr(dat,method_list = c('mr_ivw_fe',
                                'mr_egger_regression',
                                'mr_weighted_median',
                                'mr_weighted_mode'))
mr_heterogeneity(dat)
mr_pleiotropy_test(dat)
df <- df %>% bind_rows(
  data.frame(
    exposure  = 'MRI-HFC',
    outcome   = "Vascular dementia",
    method    = mr_res$method,
    nsnp      = mr_res$nsnp,
    b         = mr_res$b,
    se        = mr_res$se,
    pval      = mr_res$pval
  ))
# d <- mr_scatter_plot3(mr_res,dat)
# d2 <- d$LKFN1o.4TmsH2+labs(title = 'Hypercholesterolemia')
rm(list = c('VDE','VDE2'))

###Delirium due to conditions classified elsewhere
DEL <- fread('finngen_R7_F5_DELIRIUM.gz')
DEL2 <- format_data(DEL,
                    type = 'outcome',
                    snps = hfc_mr$snp,
                    snp_col =  'rsids',
                    beta_col = 'beta',
                    se_col = 'sebeta',
                    eaf_col = 'af_alt',
                    effect_allele_col = 'alt',
                    other_allele_col = 'ref',
                    pval_col = 'pval',
                    chr_col = '#chrom',
                    pos_col = 'pos')
dat <- harmonise_data(hfc_mr,DEL2)
write.csv(dat,'mr_Delirium.csv',row.names = F)
dat <- dat %>% mutate(mr_keep=T)
mr_res = mr(dat,method_list = c('mr_ivw_fe',
                                'mr_egger_regression',
                                'mr_weighted_median',
                                'mr_weighted_mode'))
mr_heterogeneity(dat)
mr_pleiotropy_test(dat)
df <- df %>% bind_rows(
  data.frame(
    exposure  = 'MRI-HFC',
    outcome   = "Delirium",
    method    = mr_res$method,
    nsnp      = mr_res$nsnp,
    b         = mr_res$b,
    se        = mr_res$se,
    pval      = mr_res$pval
  ))
rm(list = c('DEL','DEL2'))

###Neurological disorders
ND <- fread('finngen_R7_G6_NEURO.gz')
ND2 <- format_data(ND,
                   type = 'outcome',
                   snps = hfc_mr$snp,
                   snp_col =  'rsids',
                   beta_col = 'beta',
                   se_col = 'sebeta',
                   eaf_col = 'af_alt',
                   effect_allele_col = 'alt',
                   other_allele_col = 'ref',
                   pval_col = 'pval',
                   chr_col = '#chrom',
                   pos_col = 'pos')
dat <- harmonise_data(hfc_mr,ND2)
write.csv(dat,'mr_Neurological disorders.csv',row.names = F)
dat <- dat %>% mutate(mr_keep=T)
mr_res = mr(dat,method_list = c('mr_ivw_fe',
                                'mr_egger_regression',
                                'mr_weighted_median',
                                'mr_weighted_mode'))
mr_heterogeneity(dat)
mr_pleiotropy_test(dat)
df <- df %>% bind_rows(
  data.frame(
    exposure  = 'MRI-HFC',
    outcome   = "Neurological disorders",
    method    = mr_res$method,
    nsnp      = mr_res$nsnp,
    b         = mr_res$b,
    se        = mr_res$se,
    pval      = mr_res$pval
  ))
rm(list = c('ND','ND2'))

###phobia
PB <- fread('finngen_R7_F5_SOCPHOB.gz')
PB2 <- format_data(PB,
                   type = 'outcome',
                   snps = hfc_mr$snp,
                   snp_col =  'rsids',
                   beta_col = 'beta',
                   se_col = 'sebeta',
                   eaf_col = 'af_alt',
                   effect_allele_col = 'alt',
                   other_allele_col = 'ref',
                   pval_col = 'pval',
                   chr_col = '#chrom',
                   pos_col = 'pos')
dat <- harmonise_data(hfc_mr,PB2)
write.csv(dat,'mr_phobia.csv',row.names = F)
dat <- dat %>% mutate(mr_keep=T)
mr_res = mr(dat,method_list = c('mr_ivw_fe',
                                'mr_egger_regression',
                                'mr_weighted_median',
                                'mr_weighted_mode'))
mr_heterogeneity(dat)
mr_pleiotropy_test(dat)
df <- df %>% bind_rows(
  data.frame(
    exposure  = 'MRI-HFC',
    outcome   = "Phobia",
    method    = mr_res$method,
    nsnp      = mr_res$nsnp,
    b         = mr_res$b,
    se        = mr_res$se,
    pval      = mr_res$pval
  ))
rm(list = c('PB','PB2'))

###Alcoholic liver damage
ALD <- fread('finngen_R7_ALCOLIVER.gz')
ALD2 <- format_data(ALD,
                    type = 'outcome',
                    snps = hfc_mr$snp,
                    snp_col =  'rsids',
                    beta_col = 'beta',
                    se_col = 'sebeta',
                    eaf_col = 'af_alt',
                    effect_allele_col = 'alt',
                    other_allele_col = 'ref',
                    pval_col = 'pval',
                    chr_col = '#chrom',
                    pos_col = 'pos')
dat <- harmonise_data(hfc_mr,ALD2)
write.csv(dat,'mr_ALD.csv',row.names = F)
dat <- dat %>% mutate(mr_keep=T)
mr_res = mr(dat,method_list = c('mr_ivw_fe',
                                'mr_egger_regression',
                                'mr_weighted_median',
                                'mr_weighted_mode'))
mr_heterogeneity(dat)
mr_pleiotropy_test(dat)
df <- df %>% bind_rows(
  data.frame(
    exposure  = 'MRI-HFC',
    outcome   = "Alcoholic liver damage",
    method    = mr_res$method,
    nsnp      = mr_res$nsnp,
    b         = mr_res$b,
    se        = mr_res$se,
    pval      = mr_res$pval
  ))
e <- mr_scatter_plot3(mr_res,dat)
e2 <- e$LKFN1o.yLzq3X+labs(title = 'Alcoholic liver damage')
rm(list = c('ALD','ALD2'))

##Encephalopathy
ENC <- fread('finngen_R7_G6_ENCEPATH.gz')
ENC2 <- format_data(ENC,
                    type = 'outcome',
                    snps = hfc_mr$snp,
                    snp_col =  'rsids',
                    beta_col = 'beta',
                    se_col = 'sebeta',
                    eaf_col = 'af_alt',
                    effect_allele_col = 'alt',
                    other_allele_col = 'ref',
                    pval_col = 'pval',
                    chr_col = '#chrom',
                    pos_col = 'pos')
dat <- harmonise_data(hfc_mr,ENC2)
write.csv(dat,'mr_Encephalopathy.csv',row.names = F)
dat <- dat %>% mutate(mr_keep=T)
mr_res = mr(dat,method_list = c('mr_ivw_fe',
                                'mr_egger_regression',
                                'mr_weighted_median',
                                'mr_weighted_mode'))
mr_heterogeneity(dat)
mr_pleiotropy_test(dat)
df <- df %>% bind_rows(
  data.frame(
    exposure  = 'MRI-HFC',
    outcome   = "Encephalopathy",
    method    = mr_res$method,
    nsnp      = mr_res$nsnp,
    b         = mr_res$b,
    se        = mr_res$se,
    pval      = mr_res$pval
  ))
rm(list = c('ENC','ENC2'))
##Essential hypertension
HYP <- fread('finngen_R7_I9_HYPTENSESS.gz')
HYP2 <- format_data(HYP,
                    type = 'outcome',
                    snps = hfc_mr$snp,
                    snp_col =  'rsids',
                    beta_col = 'beta',
                    se_col = 'sebeta',
                    eaf_col = 'af_alt',
                    effect_allele_col = 'alt',
                    other_allele_col = 'ref',
                    pval_col = 'pval',
                    chr_col = '#chrom',
                    pos_col = 'pos')
dat <- harmonise_data(hfc_mr,HYP2)
write.csv(dat,'mr_EssentialHypertension.csv',row.names = F)
dat <- dat %>% mutate(mr_keep=T)
mr_res = mr(dat,method_list = c('mr_ivw_fe',
                                'mr_egger_regression',
                                'mr_weighted_median',
                                'mr_weighted_mode'))
mr_heterogeneity(dat)
mr_pleiotropy_test(dat)
df <- df %>% bind_rows(
  data.frame(
    exposure  = 'MRI-HFC',
    outcome   = "Essential hypertension",
    method    = mr_res$method,
    nsnp      = mr_res$nsnp,
    b         = mr_res$b,
    se        = mr_res$se,
    pval      = mr_res$pval
  ))
rm(list = c('HYP','HYP2'))
##CAD
cad <- fread('cad_european.tsv.gz')
cad2 <- format_data(cad,
                    type = 'outcome',
                    snps = hfc_mr$snp,
                    snp_col =  'name',
                    beta_col = 'beta',
                    se_col = 'standard_error',
                    eaf_col = 'effect_allele_frequency',
                    effect_allele_col = 'effect_allele',
                    other_allele_col = 'other_allele',
                    pval_col = 'p_value',
                    chr_col = 'chromosome',
                    pos_col = 'base_pair_location')
dat <- harmonise_data(hfc_mr,cad2)
write.csv(dat,'mr_CAD.csv',row.names = F)
dat <- dat %>% mutate(mr_keep=T)
mr_res = mr(dat,method_list = c('mr_ivw_fe',
                                'mr_egger_regression',
                                'mr_weighted_median',
                                'mr_weighted_mode'))
mr_heterogeneity(dat)
mr_pleiotropy_test(dat)
df <- df %>% bind_rows(
  data.frame(
    exposure  = 'MRI-HFC',
    outcome   = "Coronary artery disease",
    method    = mr_res$method,
    nsnp      = mr_res$nsnp,
    b         = mr_res$b,
    se        = mr_res$se,
    pval      = mr_res$pval
  ))
f <- mr_scatter_plot3(mr_res,dat)
f2 <- f$LKFN1o.Bg4hI4+labs(title = 'Coronary artery disease')
rm(list = c('cad','cad2'))

##Nonspecific chest pain
CP <- fread('finngen_R7_R18_PAIN_THROAT_CHEST.gz')
CP2 <- format_data(CP,
                   type = 'outcome',
                   snps = hfc_mr$snp,
                   snp_col =  'rsids',
                   beta_col = 'beta',
                   se_col = 'sebeta',
                   eaf_col = 'af_alt',
                   effect_allele_col = 'alt',
                   other_allele_col = 'ref',
                   pval_col = 'pval',
                   chr_col = '#chrom',
                   pos_col = 'pos')
dat <- harmonise_data(hfc_mr,CP2)
write.csv(dat,'mr_chest_pain.csv',row.names = F)
dat <- dat %>% mutate(mr_keep=T)
mr_res = mr(dat,method_list = c('mr_ivw_fe',
                                'mr_egger_regression',
                                'mr_weighted_median',
                                'mr_weighted_mode'))
mr_heterogeneity(dat)
mr_pleiotropy_test(dat)
df <- df %>% bind_rows(
  data.frame(
    exposure  = 'MRI-HFC',
    outcome   = "Nonspecific chest pain",
    method    = mr_res$method,
    nsnp      = mr_res$nsnp,
    b         = mr_res$b,
    se        = mr_res$se,
    pval      = mr_res$pval
  ))
g <- mr_scatter_plot3(mr_res,dat)
g2 <- g$LKFN1o.sXxTSj+labs(title = 'Nonspecific chest pain')
rm(list = c('CP','CP2'))

##liver cirrhosis
LCi <- fread('finngen_R7_CIRRHOSIS_BROAD.gz')
LCi2 <- format_data(LCi,
                    type = 'outcome',
                    snps = hfc_mr$snp,
                    snp_col =  'rsids',
                    beta_col = 'beta',
                    se_col = 'sebeta',
                    eaf_col = 'af_alt',
                    effect_allele_col = 'alt',
                    other_allele_col = 'ref',
                    pval_col = 'pval',
                    chr_col = '#chrom',
                    pos_col = 'pos')
dat <- harmonise_data(hfc_mr,LCi2)
write.csv(dat,'mr_cirrhosis.csv',row.names = F)
dat <- dat %>% mutate(mr_keep=T)
mr_res = mr(dat,method_list = c('mr_ivw_fe',
                                'mr_egger_regression',
                                'mr_weighted_median',
                                'mr_weighted_mode'))
mr_heterogeneity(dat)
mr_pleiotropy_test(dat)
df <- df %>% bind_rows(
  data.frame(
    exposure  = 'MRI-HFC',
    outcome   = "Cirrhosis",
    method    = mr_res$method,
    nsnp      = mr_res$nsnp,
    b         = mr_res$b,
    se        = mr_res$se,
    pval      = mr_res$pval
  ))
h <- mr_scatter_plot3(mr_res,dat)
h2 <- h$LKFN1o.6SsrXD+labs(title = 'Cirrhosis')
rm(list = c('LCi','LCi2'))

###Other disorders of stomach and duodenum
STO <- fread('finngen_R7_K11_OTHDISSTOMDUOD.gz')
STO2 <- format_data(STO,
                    type = 'outcome',
                    snps = hfc_mr$snp,
                    snp_col =  'rsids',
                    beta_col = 'beta',
                    se_col = 'sebeta',
                    eaf_col = 'af_alt',
                    effect_allele_col = 'alt',
                    other_allele_col = 'ref',
                    pval_col = 'pval',
                    chr_col = '#chrom',
                    pos_col = 'pos')
dat <- harmonise_data(hfc_mr,STO2)
write.csv(dat,'mr_Other disorders of stomach and duodenum.csv',row.names = F)
dat <- dat %>% mutate(mr_keep=T)
mr_res = mr(dat,method_list = c('mr_ivw_fe',
                                'mr_egger_regression',
                                'mr_weighted_median',
                                'mr_weighted_mode'))
mr_heterogeneity(dat)
mr_pleiotropy_test(dat)
df <- df %>% bind_rows(
  data.frame(
    exposure  = 'MRI-HFC',
    outcome   = "Other disorders of stomach and duodenum",
    method    = mr_res$method,
    nsnp      = mr_res$nsnp,
    b         = mr_res$b,
    se        = mr_res$se,
    pval      = mr_res$pval
  ))
rm(list = c('STO','STO2'))

##Ascites
ASC <- fread('finngen_R7_R18_ASCITES.gz')
ASC2 <- format_data(ASC,
                    type = 'outcome',
                    snps = hfc_mr$snp,
                    snp_col =  'rsids',
                    beta_col = 'beta',
                    se_col = 'sebeta',
                    eaf_col = 'af_alt',
                    effect_allele_col = 'alt',
                    other_allele_col = 'ref',
                    pval_col = 'pval',
                    chr_col = '#chrom',
                    pos_col = 'pos')
dat <- harmonise_data(hfc_mr,ASC2)
write.csv(dat,'mr_ascites.csv',row.names = F)
dat <- dat %>% mutate(mr_keep=T)
mr_res = mr(dat,method_list = c('mr_ivw_fe',
                                'mr_egger_regression',
                                'mr_weighted_median',
                                'mr_weighted_mode'))
mr_heterogeneity(dat)
mr_pleiotropy_test(dat)
df <- df %>% bind_rows(
  data.frame(
    exposure  = 'MRI-HFC',
    outcome   = "Ascites",
    method    = mr_res$method,
    nsnp      = mr_res$nsnp,
    b         = mr_res$b,
    se        = mr_res$se,
    pval      = mr_res$pval
  ))
k <- mr_scatter_plot3(mr_res,dat)
k2 <- k$LKFN1o.rFc2vH+labs(title = 'Ascites')
rm(list = c('ASC','ASC2'))


##Splenomegaly
SP <- fread('finngen_R7_R18_HEPATOM_SPLENO_NOT_ELSEW_CLASSIFIED.gz')
SP2 <- format_data(SP,
                   type = 'outcome',
                   snps = hfc_mr$snp,
                   snp_col =  'rsids',
                   beta_col = 'beta',
                   se_col = 'sebeta',
                   eaf_col = 'af_alt',
                   effect_allele_col = 'alt',
                   other_allele_col = 'ref',
                   pval_col = 'pval',
                   chr_col = '#chrom',
                   pos_col = 'pos')
dat <- harmonise_data(hfc_mr,SP2)
write.csv(dat,'mr_Splenomegaly.csv',row.names = F)
dat <- dat %>% mutate(mr_keep=T)
mr_res = mr(dat,method_list = c('mr_ivw_fe',
                                'mr_egger_regression',
                                'mr_weighted_median',
                                'mr_weighted_mode'))
mr_heterogeneity(dat)
mr_pleiotropy_test(dat)
df <- df %>% bind_rows(
  data.frame(
    exposure  = 'MRI-HFC',
    outcome   = "Splenomegaly",
    method    = mr_res$method,
    nsnp      = mr_res$nsnp,
    b         = mr_res$b,
    se        = mr_res$se,
    pval      = mr_res$pval
  ))
rm(list = c('SP','SP2'))
write.csv(df,'mr_mrihfc_results.csv',row.names = F)

####plot for MR results
plot_grid(a2,b2,cc2,d2,e2,f2,g2,h2,k2,ncol = 3)
mrs <- df
mrs <- mrs %>% mutate(
  OR = exp(b),
  lwr= exp(b-1.96*se),
  lwr= exp(b+1.96*se),
  fdr=p.adjust(pval,'fdr')
)
mrs <- mrs %>% mutate(
  outcome =  factor(outcome,levels = unique(mrs$outcome)),
  or_group =  cut(OR,breaks =  c(0,0.6,0.8,1,1.5,2,2.5,3,100),
                  labels = c('<0.60','0.60-0.80','0.81-1.00','1.01-1.50',
                             '1.51-2.00','2.01-2.50','2.51-3.00','>3.00'),
                  include_lowest = T,right = T)
)

mrs <- mrs %>% mutate(
  sigs = ifelse(fdr<0.001,'***',ifelse(fdr<0.01,'**',ifelse(fdr<0.05,'*',''))),
  method =  factor(method,levels = unique(mrs$method))
)

a <- ggplot(mrs,aes(method,outcome))+
  geom_tile(aes(fill = or_group),color='gray20')+
  scale_fill_manual(name='Odds ratio',
                    values = c('#08519c','#6baed6','#c6dbef',
                               '#fee0d2','#fcbba1','#fc9272','#ef3b2c','#cb181d'))+
  geom_text(aes(label = sprintf("%0.2f",OR)))+
  geom_text(aes(label = sigs),vjust=1.8)+
  geom_vline(xintercept = c(1.5,2.5,3.5,4.5),color='gray20')+
  geom_hline(yintercept = seq(0.5,19.5,by=1),color='gray20')+
  theme_test(base_size = 12)+labs(x=NULL,y=NULL)+
  scale_y_discrete(expand = c(0,0),limits=rev)+
  scale_x_discrete(expand = c(0,0),labels= c('IVW','MR Egger','Weighted median','Weighted mode'))+
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        axis.text = element_text(size=12))

h_test <- readxl::read_excel('GWAS source.xlsx',sheet = 2)
h_test <- h_test %>% na.omit()
h_test <- h_test %>% mutate(
  fdr_heterogeneitye = p.adjust(`P-Q`,'fdr'),
  fdr_pleiotropy     = p.adjust(`Egger-P`,'fdr')
)

h_test <- h_test %>% pivot_longer(cols =  starts_with("fdr"),
                                  names_to = 'test',
                                  values_to = 'FDR')
h_test <- h_test %>% mutate(
  phenotype = factor(phenotype,levels= unique(h_test$phenotype)),
  sigs = ifelse(FDR<0.001,'***',ifelse(FDR<0.01,'**',ifelse(FDR<0.05,'*','')))
)

b <- ggplot(h_test,aes(test,phenotype))+
  geom_tile(aes(fill = sigs),color='gray20',fill='white')+
  geom_text(aes(label = scales::scientific(FDR,digits = 2)),size=3)+
  geom_text(aes(label = sigs),vjust=1.8)+
  geom_vline(xintercept = c(1.5,2.5,3.5,4.5),color='gray20')+
  geom_hline(yintercept = seq(0.5,19.5,by=1),color='gray20')+
  theme_test(base_size = 12)+labs(x=NULL,y=NULL)+
  scale_y_discrete(expand = c(0,0),limits=rev)+
  scale_x_discrete(expand = c(0,0),labels= c('Heterogeneity','Horizontal\npleiotropy'))+
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        axis.text = element_text(size=12))
plot_grid(a+theme(legend.position = 'none'),
          b,ncol = 2,align = 'h',
          rel_widths = c(1,.8))

###---MVMR--------------------------------------------------------------------------
library(MVMR)
library(TwoSampleMR)
####finding the significant SNPs for all GWAS
bmi <- fread('Meta-analysis_Locke_et_al+UKBiobank_2018_UPDATED_BMI.txt.gz')
smoking <- fread('D:\\GWAS of smoking and alcohol\\CigarettesPerDay.txt.gz')
alcohol <- fread('D:\\GWAS of smoking and alcohol\\DrinksPerWeek.txt.gz')
hfc <- read_rds('gwas0_PDFF_mhtest0426.rds')
whr2 <- extract_instruments(outcomes = 'ieu-a-72',clump = F)
whr2 <- whr2 %>% select(SNP)
bmi2 <- bmi %>% filter(P<5e-8) %>% select(SNP)
smoking2 <- smoking %>% filter(PVALUE<5e-8) %>% select(RSID) %>% rename(SNP=RSID)
alcohol2 <- alcohol %>% filter(PVALUE<5e-8) %>% select(RSID)%>% rename(SNP=RSID)
hfc2 <- hfc %>% filter(pval<5e-8) %>% select(snp)%>% rename(SNP=snp)
union_snps <- hfc2 %>% union(bmi2) %>% union(smoking2) %>% union(alcohol2) %>% union(whr2)
union_snps <- union_snps %>% na.omit()
####clumping the union snps using one of GWAS
hfc2c <- hfc %>% filter(snp %in% union_snps$SNP) %>% select(chr,pos,pval,snp) 
hfc2c <- hfc2c %>% rename(chr_name=chr,chrom_start=pos,pval.exposure=pval,SNP=snp) %>%
  clump_data(clump_r2 = 0.01)

###extracting data and combining them together
bmi3 <- bmi %>% filter(SNP %in% hfc2c$SNP) %>% select(SNP,BETA,SE) %>% rename(beta.bmi=BETA,se.bmi=SE)
smoking3 <- smoking %>% filter(RSID %in% hfc2c$SNP)%>% select(RSID,BETA,SE) %>% 
  rename(SNP = RSID,beta.smoking=BETA,se.smoking=SE)
alcohol3 <- alcohol %>% filter(RSID %in% hfc2c$SNP)%>% select(RSID,BETA,SE) %>% 
  rename(SNP = RSID,beta.alcohol=BETA,se.alcohol=SE)
hfc3 <- hfc %>% filter(snp %in% hfc2c$SNP) %>% select(beta,se,snp) %>% 
  rename(SNP=snp, beta.hfc=beta,se.hfc=se)
whr3 <- extract_outcome_data(snps = hfc2c$SNP,outcomes = 'ieu-a-72')
whr3 <- whr3 %>% select(SNP,beta.outcome,se.outcome) %>% rename(beta.whr=beta.outcome,se.whr=se.outcome)

exposure.data <- bmi3 %>% left_join(smoking3) %>% left_join(alcohol3) %>% left_join(hfc3) %>% left_join(whr3) %>%
  na.omit()


LC <- fread('finngen_R7_C3_LIVER_INTRAHEPATIC_BILE_DUCTS_EXALLC.gz')
LC2 <- LC %>% filter(rsids %in% exposure.data$SNP) %>% select(beta,sebeta,rsids) %>%
  rename(SNP=rsids,beta.out=beta,se.out=sebeta)

rawdat_mvmr <- exposure.data %>% left_join(LC2) %>% na.omit()
F.data <- format_mvmr(BXGs = rawdat_mvmr[,c(2,4,6,8,10)],
                      BYG = rawdat_mvmr[,12],
                      seBXGs = rawdat_mvmr[,c(3,5,7,9,11)],
                      seBYG = rawdat_mvmr[,13],
                      RSID = rawdat_mvmr[,1])
sres <- strength_mvmr(r_input = F.data, gencov = 0)
pres <- pleiotropy_mvmr(r_input = F.data, gencov = 0)
res <- ivw_mvmr(r_input = F.data)
df <- as.data.frame(res) %>% mutate(
  exposure = c('BMI','Smoking','Alcohol','HFC','WHR'),
  outcome  = 'Primary liver cancer'
)

CHP <- fread('finngen_R7_K11_CHRONHEP.gz')
CHR2 <- CHP %>% filter(rsids %in% exposure.data$SNP) %>% select(beta,sebeta,rsids) %>%
  rename(SNP=rsids,beta.CHP=beta,se.CHP=sebeta)

rawdat_mvmr <- exposure.data %>% left_join(CHR2) %>% na.omit()
F.data <- format_mvmr(BXGs = rawdat_mvmr[,c(2,4,6,8,10)],
                      BYG = rawdat_mvmr[,12],
                      seBXGs = rawdat_mvmr[,c(3,5,7,9,11)],
                      seBYG = rawdat_mvmr[,13],
                      RSID = rawdat_mvmr[,1])
sres <- strength_mvmr(r_input = F.data, gencov = 0)
pres <- pleiotropy_mvmr(r_input = F.data, gencov = 0)
res <- ivw_mvmr(r_input = F.data)
df <- df %>% bind_rows(as.data.frame(res) %>% mutate(
  exposure = c('BMI','Smoking','Alcohol','HFC','WHR'),
  outcome  = 'Hepatitis NOS'
))

t2d <- fread('Mahajan.NatGenet2018b.T2D-noUKBB.European.txt')
hfc2c <- hfc2c %>% mutate(SNP2 = paste0(chr_name,":",chrom_start))
t2d2 <- t2d %>% filter(SNP %in% hfc2c$SNP2) %>% select(SNP,Beta,SE)
t2d2 <- t2d2 %>% rename(SNP2=SNP) %>% left_join(hfc2c %>% select(SNP2,SNP),by = c('SNP2'))
t2d2 <- t2d2 %>% select(-1) %>% rename(
  beta.out=Beta,se.out=SE
)
rawdat_mvmr <- exposure.data %>% left_join(t2d2) %>% na.omit()
F.data <- format_mvmr(BXGs = rawdat_mvmr[,c(2,4,6,8,10)],
                      BYG = rawdat_mvmr[,12],
                      seBXGs = rawdat_mvmr[,c(3,5,7,9,11)],
                      seBYG = rawdat_mvmr[,13],
                      RSID = rawdat_mvmr[,1])
sres <- strength_mvmr(r_input = F.data, gencov = 0)
pres <- pleiotropy_mvmr(r_input = F.data, gencov = 0)
res <- ivw_mvmr(r_input = F.data)
df <- df %>% bind_rows(as.data.frame(res) %>% mutate(
  exposure = c('BMI','Smoking','Alcohol','HFC','WHR'),
  outcome  = 'Type 2 diabetes'
))

HLA <- fread('DYSLIPID_GERA.txt.gz')
HLA2 <- HLA %>% filter(name %in% exposure.data$SNP) %>% select(frequentist_add_beta_1,frequentist_add_se_1,name) %>%
  rename(SNP=name,beta.out=frequentist_add_beta_1,se.out=frequentist_add_se_1)

rawdat_mvmr <- exposure.data %>% left_join(HLA2) %>% na.omit()
F.data <- format_mvmr(BXGs = rawdat_mvmr[,c(2,4,6,8,10)],
                      BYG = rawdat_mvmr[,12],
                      seBXGs = rawdat_mvmr[,c(3,5,7,9,11)],
                      seBYG = rawdat_mvmr[,13],
                      RSID = rawdat_mvmr[,1])
sres <- strength_mvmr(r_input = F.data, gencov = 0)
pres <- pleiotropy_mvmr(r_input = F.data, gencov = 0)
res <- ivw_mvmr(r_input = F.data)
df <- df %>% bind_rows(as.data.frame(res) %>% mutate(
  exposure = c('BMI','Smoking','Alcohol','HFC','WHR'),
  outcome  = 'Hypercholesterolemia'
))

ALD <- fread('finngen_R7_ALCOLIVER.gz')
ALD2 <- ALD %>% filter(rsids %in% exposure.data$SNP) %>% select(beta,sebeta,rsids) %>%
  rename(SNP=rsids,beta.out=beta,se.out=sebeta)

rawdat_mvmr <- exposure.data %>% left_join(ALD2) %>% na.omit()
F.data <- format_mvmr(BXGs = rawdat_mvmr[,c(2,4,6,8,10)],
                      BYG = rawdat_mvmr[,12],
                      seBXGs = rawdat_mvmr[,c(3,5,7,9,11)],
                      seBYG = rawdat_mvmr[,13],
                      RSID = rawdat_mvmr[,1])
sres <- strength_mvmr(r_input = F.data, gencov = 0)
pres <- pleiotropy_mvmr(r_input = F.data, gencov = 0)
res <- ivw_mvmr(r_input = F.data)
df <- df %>% bind_rows(as.data.frame(res) %>% mutate(
  exposure = c('BMI','Smoking','Alcohol','HFC','WHR'),
  outcome  = 'Alcoholic liver damage'
))

cad <- fread('cad_european.tsv.gz')
cad2 <- cad %>% filter(name %in% exposure.data$SNP) %>% select(beta,standard_error,name) %>%
  rename(SNP=name,beta.out=beta,se.out=standard_error)

rawdat_mvmr <- exposure.data %>% left_join(cad2) %>% na.omit()
F.data <- format_mvmr(BXGs = rawdat_mvmr[,c(2,4,6,8,10)],
                      BYG = rawdat_mvmr[,12],
                      seBXGs = rawdat_mvmr[,c(3,5,7,9,11)],
                      seBYG = rawdat_mvmr[,13],
                      RSID = rawdat_mvmr[,1])
sres <- strength_mvmr(r_input = F.data, gencov = 0)
pres <- pleiotropy_mvmr(r_input = F.data, gencov = 0)
res <- ivw_mvmr(r_input = F.data)
df <- df %>% bind_rows(as.data.frame(res) %>% mutate(
  exposure = c('BMI','Smoking','Alcohol','HFC','WHR'),
  outcome  = 'Coronary artery disease'
))

CP <- fread('finngen_R7_R18_PAIN_THROAT_CHEST.gz')
CP2 <- CP %>% filter(rsids %in% exposure.data$SNP) %>% select(beta,sebeta,rsids) %>%
  rename(SNP=rsids,beta.out=beta,se.out=sebeta)

rawdat_mvmr <- exposure.data %>% left_join(CP2) %>% na.omit()
F.data <- format_mvmr(BXGs = rawdat_mvmr[,c(2,4,6,8,10)],
                      BYG = rawdat_mvmr[,12],
                      seBXGs = rawdat_mvmr[,c(3,5,7,9,11)],
                      seBYG = rawdat_mvmr[,13],
                      RSID = rawdat_mvmr[,1])
sres <- strength_mvmr(r_input = F.data, gencov = 0)
pres <- pleiotropy_mvmr(r_input = F.data, gencov = 0)
res <- ivw_mvmr(r_input = F.data)
df <- df %>% bind_rows(as.data.frame(res) %>% mutate(
  exposure = c('BMI','Smoking','Alcohol','HFC','WHR'),
  outcome  = 'Nonspecific chest pain'
))

LCi <- fread('finngen_R7_CIRRHOSIS_BROAD.gz')
LCi2 <- LCi %>% filter(rsids %in% exposure.data$SNP) %>% select(beta,sebeta,rsids) %>%
  rename(SNP=rsids,beta.out=beta,se.out=sebeta)

rawdat_mvmr <- exposure.data %>% left_join(LCi2) %>% na.omit()
F.data <- format_mvmr(BXGs = rawdat_mvmr[,c(2,4,6,8,10)],
                      BYG = rawdat_mvmr[,12],
                      seBXGs = rawdat_mvmr[,c(3,5,7,9,11)],
                      seBYG = rawdat_mvmr[,13],
                      RSID = rawdat_mvmr[,1])
sres <- strength_mvmr(r_input = F.data, gencov = 0)
pres <- pleiotropy_mvmr(r_input = F.data, gencov = 0)
res <- ivw_mvmr(r_input = F.data)
df <- df %>% bind_rows(as.data.frame(res) %>% mutate(
  exposure = c('BMI','Smoking','Alcohol','HFC','WHR'),
  outcome  = 'Liver cirrhosis'
))

ASC <- fread('finngen_R7_R18_ASCITES.gz')
ASC2 <- ASC %>% filter(rsids %in% exposure.data$SNP) %>% select(beta,sebeta,rsids) %>%
  rename(SNP=rsids,beta.out=beta,se.out=sebeta)

rawdat_mvmr <- exposure.data %>% left_join(ASC2) %>% na.omit()
F.data <- format_mvmr(BXGs = rawdat_mvmr[,c(2,4,6,8,10)],
                      BYG = rawdat_mvmr[,12],
                      seBXGs = rawdat_mvmr[,c(3,5,7,9,11)],
                      seBYG = rawdat_mvmr[,13],
                      RSID = rawdat_mvmr[,1])
sres <- strength_mvmr(r_input = F.data, gencov = 0)
pres <- pleiotropy_mvmr(r_input = F.data, gencov = 0)
res <- ivw_mvmr(r_input = F.data)
df <- df %>% bind_rows(as.data.frame(res) %>% mutate(
  exposure = c('BMI','Smoking','Alcohol','HFC','WHR'),
  outcome  = 'Ascites'
))

write.csv(df,'MVMR_MRI_HFC.csv',row.names = F)

mvmr_res <- read.csv('MVMR_MRI_HFC.csv')
mvmr_res <- mvmr_res %>% mutate(
  exposure = factor(exposure,levels = rev(c('HFC','BMI','Smoking','Alcohol','WHR')),
                    labels = rev(c('Hepatic fat content','Body mass index','Smoking',
                               'Alcohol drinking','Waist-hip-ratio'))),
  outcome  = factor(outcome,levels = rev(c('Hepatitis NOS','Primary liver cancer','Type 2 diabetes',
                                           'Hypercholesterolemia','Alcoholic liver damage',
                                           'Coronary artery disease','Nonspecific chest pain',
                                           'Liver cirrhosis','Ascites'))),
  OR      = exp(Estimate),
  lwr     = exp(Estimate-1.96*Std..Error),
  upr     = exp(Estimate+1.96*Std..Error),
  p       = Pr...t..
)

ggplot(mvmr_res,aes(OR,outcome,color = exposure))+
  geom_vline(xintercept = 1,linetype = 2)+
  geom_hline(yintercept = seq(1.5,9.5,by=1),color='gray70',linetype=2)+
  geom_errorbarh(aes(xmin=lwr,xmax=upr),position = position_dodge(width=.8),height=0)+
  geom_point(size=3,position = position_dodge(width=.8))+
  paletteer::scale_colour_paletteer_d("ggthemes::few_Dark",
                                      breaks  = c('Hepatic fat content','Body mass index','Smoking',
                                                  'Alcohol drinking','Waist-hip-ratio'))+
  scale_x_continuous(limits = c(0,5),breaks = seq(.5,3.5,by=.5))+
  scale_y_discrete(expand = expansion(add = c(0.5,0.75)))+
  # geom_text(aes(x = 4,y=outcome,label = `OR (95% CI)`,group=exposure),
  #           position=position_dodge(width=.8),
  #           color='black')+
  # annotate('text',x = rep(3.,4), y= c(4,3,2,1),label = c(877,919,887,930),
  #          color='black')+
  geom_text(aes(x = 4.6,y=outcome,label = scales::scientific(Pr...t..,digits = 2),
                group=exposure),color='black',position=position_dodge(width=.8))+
  # annotate('text',x = 4,y=4.6,label = 'OR (95% CI)',color='black',fontface='bold')+
  # annotate('text',x = 3.,y=4.6,label = 'Valid SNP',color='black',fontface='bold')+
  annotate('text',x = 4.6,y=9.65,label = 'P value',color='black',fontface='bold')+
  labs(x = 'Odds ratio (95% CI)',y=NULL)+
  theme_bw()+
  guides(color = guide_legend(ncol = 2))+
  theme(axis.text = element_text(size=10,color='black'),
        legend.key = element_rect(fill = NA),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        legend.position = 'top')


####---mediation analysis--------------------------------------------------------------
####run on LINUX
library(mediation)
tempData2 <- fread('data_for_mediation.csv')
disease_list <- list(c(530.2,571.81,571.8,571.5,571.51),
                     c(155,155.1),
                     c(70.9),
                     c(250.2),
                     c(272.11),
                     c(317.11),
                     c(411.3,411.2,411.8),
                     c(418),
                     c(572))
outcomes <- c('Liver cirrhosis','PLC','Hepatitis NOS','T2D','Hypercholesterolemia',
              'ALD','CAD','Chest pain','Ascites')
df2 <- data.frame()
for(i in 1:9){
  case <- tempData2 %>% filter(PheCode %in% disease_list[[i]]) %>% 
    distinct(eid,.keep_all = TRUE) %>%mutate(pheno =1)
  phenotypeGroup <- tempData2 %>% filter(PheCode== disease_list[[i]][1]) %>% select(PhenotypeGroup)
  control <- tempData2 %>% filter(PhenotypeGroup != phenotypeGroup$PhenotypeGroup[1]) %>% 
    filter(!eid %in% case$eid) %>%
    distinct(eid,.keep_all = TRUE) %>% 
    mutate(pheno=0)
  dat1 <- case %>% bind_rows(control)
  bios <- names(tempData2)[15:81]
  if(i==4){
    dat1 <- dat1 %>% filter(loweringINS==0)
    for(j in 1:67){
      biomarker = get(bios[j],pos = dat1)
      med.fit <- lm(biomarker~HFCGRS+sex+age+smokingStatus+averageTotalIncome+aolIntakeFreq+
                      no.daysModeAct+education+BMI,data = dat1)
      out.fit <- glm(pheno~HFCGRS+sex+age+smokingStatus+averageTotalIncome+aolIntakeFreq+
                     no.daysModeAct+education+BMI+biomarker,data = dat1,family = binomial)
      mediating <- mediate(med.fit,out.fit,treat = 'HFCGRS',mediator = "biomarker")
      df <- data.frame(
        outcome  = outcomes[i],
        mediator = bios[j],
        pm       = mediating$n.avg,
        pm_lwr   = mediating$n.avg.ci[1],
        pm_upr   = mediating$n.avg.ci[2],
        p        = mediating$n.avg.p
      )
      df2 <- bind_rows(df2,df)
      print(paste(outcomes[i],":",j,'th biomarker was finished!!!'))
    }
  }else if(i==7){
    dat1 <- dat1 %>% filter(loweringCHL==0)
    for(j in 1:67){
      biomarker = get(bios[j],pos = dat1)
      med.fit <- lm(biomarker~HFCGRS+sex+age+smokingStatus+averageTotalIncome+aolIntakeFreq+
                      no.daysModeAct+education+BMI,data = dat1)
      out.fit <- glm(pheno~HFCGRS+sex+age+smokingStatus+averageTotalIncome+aolIntakeFreq+
                       no.daysModeAct+education+BMI+biomarker,data = dat1,family = binomial)
      mediating <- mediate(med.fit,out.fit,treat = 'HFCGRS',mediator = "biomarker")
      df <- data.frame(
        outcome  = outcomes[i],
        mediator = bios[j],
        pm       = mediating$n.avg,
        pm_lwr   = mediating$n.avg.ci[1],
        pm_upr   = mediating$n.avg.ci[2],
        p        = mediating$n.avg.p
      )
      df2 <- bind_rows(df2,df)
      print(paste(outcomes[i],":",j,'th biomarker was finished!!!'))
    }
  }else{
    for(j in 1:67){
      biomarker = get(bios[j],pos = dat1)
      med.fit <- lm(biomarker~HFCGRS+sex+age+smokingStatus+averageTotalIncome+aolIntakeFreq+
                      no.daysModeAct+education+BMI,data = dat1)
      out.fit <- glm(pheno~HFCGRS+sex+age+smokingStatus+averageTotalIncome+aolIntakeFreq+
                       no.daysModeAct+education+BMI+biomarker,data = dat1,family = binomial)
      mediating <- mediate(med.fit,out.fit,treat = 'HFCGRS',mediator = "biomarker")
      df <- data.frame(
        outcome  = outcomes[i],
        mediator = bios[j],
        pm       = mediating$n.avg,
        pm_lwr   = mediating$n.avg.ci[1],
        pm_upr   = mediating$n.avg.ci[2],
        p        = mediating$n.avg.p
      )
      df2 <- bind_rows(df2,df)
      print(paste(outcomes[i],":",j,'th biomarker was finished!!!'))
    }
  }
  write.csv(df2,'mediating.csv',row.names = F)
}

###RUN ON WINDOWS
mediating <- read.csv('mediating.csv')

###determine the significant association between HFC and mediators
hfc_to_bio <- read.csv('biomarkers_HFCGRS.csv')
hfc_to_bio <- hfc_to_bio %>% rename(hfc_to_mediator_fdr = fdr,
                                    hfc_to_mediator_beta=beta,
                                    hfc_to_mediator_se  =se)
mediating <- mediating %>% left_join(hfc_to_bio %>% 
                                       select(Biomarker,BiomarkerNames,hfc_to_mediator_fdr,
                                              hfc_to_mediator_beta, hfc_to_mediator_se),
                                     by = c('mediator'='Biomarker'))

###determine the significant association between mediators and outcomes
bios_to_out <- read.csv('biomarker-to-diseases.csv')
bios_to_out <- bios_to_out %>% rename(mediator_to_outcome_fdr = fdr,
                                      mediator_to_outcome_beta= beta,
                                      mediator_to_outcome_se  = se)

mediating <- mediating %>% left_join(bios_to_out %>% 
                                       select(biomarker,mediator_to_outcome_fdr,
                                              mediator_to_outcome_beta, 
                                              mediator_to_outcome_se,
                                              outcome),
                                     by = c('mediator'='biomarker',
                                            'outcome' ='outcome'))
mediating2 <- mediating %>% filter(
  hfc_to_mediator_fdr<0.05,
  mediator_to_outcome_fdr<0.05
)

mediating2 <- mediating2 %>% group_by(outcome) %>%
  mutate(pm.fdr = p.adjust(p,'fdr')) %>% ungroup()

mediating3 <- mediating2 %>% filter(
  pm.fdr<0.05
) %>% rename(
  pm.p=p
)

hfc_to_dis <- read.csv('HFC-to-9diseases.csv')
hfc_to_dis <- hfc_to_dis %>% rename(
  hfc_to_outcome_beta = beta,
  hfc_to_outcome_se   = se,
  hfc_to_outcome_fdr  = fdr
) %>% select(-p)

mediating3 <- mediating3 %>% left_join(hfc_to_dis)
mediating3 <- mediating3 %>% select(15,2,7,1,16:18,9,10,8,12,13,11,3,4,5,6,14)
write.csv(mediating3,'mediating_res.csv',row.names = F)

mediating3 <- mediating3 %>% mutate(
  outcome = factor(outcome,levels = unique(mediating3$outcome)[c(3,2,4,5,6,7,8,1,9)],
                   labels = c('Hepatitis NOS','Primary liver cancer','Type 2 diabetes',
                              'Hypercholesterolemia','Alcoholic liver damage',
                              'Coronary artery disease','Chest pain','Liver cirrhosis','Ascites')),
  BiomarkerNames = factor(BiomarkerNames,levels = unique(mediating3$BiomarkerNames)[c(1,40,2:15,17:22,
                                                                                      41,24,25,42,26,23,
                                                                                      27:39,16)])
)
mediating3 <- mediating3 %>% mutate(
  pmGroup = cut(pm,breaks = c(-30,-20,-10,-5,0,5,10,15,20,30,60)/100,
                labels = c('-30~-20%','-20~-10%','-10~-5%','-5~0%',
                           '0-5%','5-10%','10-15%','15-20%','20-30%','>30%'),
                include.lowest = T,right = T)
)
mediating3 <- mediating3 %>% mutate(
  twenty = ifelse(abs(pm)>=0.10,'white','black')
)

ggplot(mediating3,aes(outcome,BiomarkerNames))+
  geom_tile(aes(fill=pmGroup),color='gray20')+
  scale_fill_manual(name = 'Proportion mediated',
                    values = c('#253473FF','#08519c','#6baed6','#c6dbef',
                                                   '#fee0d2','#fcbba1','#fc9272','#ef3b2c','#cb181d','#941B0CFF'))+
  
  geom_vline(xintercept = seq(1.5,9.5,by=1),color='gray20')+
  geom_hline(yintercept = seq(0.5,42.5,by=1),color='gray20')+
  geom_text(aes(label=sprintf("%0.1f",pm*100),color=factor(twenty)))+
  scale_color_manual(values = c('black','white'),guide='none')+
  theme_test(base_size = 12)+labs(x=NULL,y=NULL)+
  scale_y_discrete(expand = c(0,0),limits=rev)+
  scale_x_discrete(expand = c(0,0))+
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        axis.text = element_text(size=11,color='black'),
        legend.key.height = unit(.1,'cm'),
        legend.spacing.y = unit(.3,'cm'))

