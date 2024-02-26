#load package
library(RODBC)

# open odbc connections - via sql authentication
conn_D1_SmartLab <- "driver={SQL Server};server=SrvSqlD1;database=D1_SmartLab;uid=smartlab_report;pwd=R3port1ng;"
conn_D1_Security <- "driver={SQL Server};server=SrvSqlD1;database=D1_Security;trusted_connection=yes;"


# Get lab.id from lab.code
GetLabID <- function(lab){
  channel <- odbcDriverConnect(connection = conn_D1_Security) 
  ds <- sqlQuery(channel, paste("SELECT TOP 1 [ID] FROM [D1_Security].[dbo].[Team] WHERE CODE = '", lab, "'", sep=""), as.is=TRUE)
  odbcClose(channel)
  return(as.numeric(ds[1,1]))}

# Get project.id from project.code
GetProjectID <- function(lab, project){
  channel <- odbcDriverConnect(connection = conn_D1_SmartLab) 
  ds <- sqlQuery(channel, paste("SELECT TOP 1 [ID] FROM [D1_SmartLab].[pick].[Project] WHERE  LabID = ", GetLabID(lab), " AND code = '", project, "' AND GcRecord = 0", sep=""), as.is=TRUE)
  odbcClose(channel)
  return(ds[1,1])}


# Load service, sampletset, analysis & sample: base properties
loadDataBaseProp <- function(lab, project, year){
  channel <- odbcDriverConnect(connection = conn_D1_SmartLab) 
  ds <- sqlQuery(channel, paste("EXEC [smart].[GetOutputBase] @Lab = ", lab, " ,@Project =  '", project, "' ,@Year = ", year, sep=""), as.is=TRUE)
  odbcClose(channel)
  ds[,"Outcome_Result"] = as.numeric(ds[,"Outcome_Result"])
  return(ds)}


# Load service: dynamic properties
loadServiceDynProp <- function(lab, project, year){
  channel <- odbcDriverConnect(connection = conn_D1_SmartLab) 
  ds <- sqlQuery(channel, paste("EXEC [smart].[GetServiceDynamicProperty] @LabID = ", GetLabID(lab), " ,@ProjectID =  '", GetProjectID(lab, project), "' ,@Year = ", year, "", sep=""), as.is=TRUE)
  odbcClose(channel)
  colnames(ds)[colnames(ds)=="RecordID"] <- "Service_ID"
  colnames(ds)[colnames(ds)=="Trip dep. date"] <- "Service_TripDepartureDate"
  colnames(ds)[colnames(ds)=="Observer"] <- "Service_Observer"
  colnames(ds)[colnames(ds)=="Vessel"] <- "Service_Vessel"
  return(ds)}


# Load sampleset: dynamic properties
loadSamplesetDynProp <- function(lab, project, year){
  channel <- odbcDriverConnect(connection = conn_D1_SmartLab) 
  ds <- sqlQuery(channel, paste("EXEC [smart].[GetSamplesetDynamicProperty] @LabID = ", GetLabID(lab), " ,@ProjectID =  '", GetProjectID(lab, project), "' ,@Year = ", year, "", sep=""), as.is=TRUE)
  odbcClose(channel)
  colnames(ds)[colnames(ds)=="RecordID"] <- "SampleSet_ID"
  colnames(ds)[colnames(ds)=="Type"] <- "SampleSet_Type"
  colnames(ds)[colnames(ds)=="Ices Division"] <- "SampleSet_IcesDivision"
  colnames(ds)[colnames(ds)=="Statistical Rectangle"] <- "SampleSet_StatisticalRectangle"  
  colnames(ds)[colnames(ds)=="Species"] <- "SampleSet_Species"
  colnames(ds)[colnames(ds)=="Fate category"] <- "SampleSet_FateCategory"
  colnames(ds)[colnames(ds)=="Folder pictures"] <- "SampleSet_FolderPictures" 
  return(ds)}


# Load sampleset_sample: dynamic properties
loadSamplesetSampleDynProp <- function(lab, project, year){
  channel <- odbcDriverConnect(connection = conn_D1_SmartLab) 
  ds <- sqlQuery(channel, paste("EXEC [smart].[GetSampleSetSampleDynamicProperty] @LabID = ", GetLabID(lab), " ,@ProjectID =  '", GetProjectID(lab, project), "' ,@Year = ", year, "", sep=""), as.is=TRUE)
  odbcClose(channel)
  colnames(ds)[colnames(ds)=="RecordID"] <- "Sample_ID"
  colnames(ds)[colnames(ds)=="Catch date"] <- "Sample_CatchDate"  
  colnames(ds)[colnames(ds)=="Length"] <- "Sample_Length"  
  colnames(ds)[colnames(ds)=="Weight"] <- "Sample_Weight"
  colnames(ds)[colnames(ds)=="Maturity"] <- "Sample_Maturity"  
  colnames(ds)[colnames(ds)=="Sex"] <- "Sample_Sex" 
  colnames(ds)[colnames(ds)=="Ices Division"] <- "Sample_IcesDivision"
  colnames(ds)[colnames(ds)=="Statistical Rectangle"] <- "Sample_StatisticalRectangle"  
  colnames(ds)[colnames(ds)=="Species"] <- "Sample_Species"
  ds[,"Sample_Length"] = as.numeric(gsub(",", ".", gsub("\\.", "", ds[,"Sample_Length"]))) # replace "," with "." before converting char into numeric
  ds[,"Sample_Weight"] = as.numeric(gsub(",", ".", gsub("\\.", "", ds[,"Sample_Weight"]))) # replace "," with "." before converting char into numeric
  return(ds)}




# Load analysis_sample: dynamic properties
loadAnalysisSampleDynProp <- function(lab, project, year){
  channel <- odbcDriverConnect(connection = conn_D1_SmartLab) 
  ds <- sqlQuery(channel, paste("EXEC [smart].[GetAnalysisSampleDynamicProperty] @LabID = ", GetLabID(lab), " ,@ProjectID =  '", GetProjectID(lab, project), "' ,@Year = ", year, "", sep=""), as.is=TRUE)
  odbcClose(channel)
  colnames(ds)[colnames(ds)=="RecordID"] <- "Sample_ID"
  colnames(ds)[colnames(ds)=="Resin No"] <- "Sample_ResinNumber"  
  colnames(ds)[colnames(ds)=="Colour No"] <- "Sample_ColourNumber"  
  colnames(ds)[colnames(ds)=="Magnification"] <- "Sample_Magnification"
  return(ds)}


# Join different datasets into 1 single dataset 
loadDataOtl <- function(lab, project, year){
  BaseProp <- loadDataBaseProp(lab, project, year)
  SE_DynProp <- loadServiceDynProp(lab, project, year)
  SS_DynProp <- loadSamplesetDynProp(lab, project, year)
  SS_S_DynProp <- loadSamplesetSampleDynProp(lab, project, year)
  A_S_DynProp <- loadAnalysisSampleDynProp(lab, project, year)
  ds <- left_join(BaseProp,SE_DynProp, by="Service_ID")
  ds <- left_join(ds,SS_DynProp, by="SampleSet_ID")
  ds <- left_join(ds,SS_S_DynProp, by="Sample_ID")
  ds <- left_join(ds,A_S_DynProp, by="Sample_ID")
  return(ds)}

# Load age reading detail: line & dot properties
loadAgeReadingDetail <- function(lab, project, year){
  channel <- odbcDriverConnect(connection = conn_D1_SmartLab) 
  ds <- sqlQuery(channel, paste("EXEC [smart].[GetAgeReadingDetail] @Lab = ", lab, " ,@Project =  '", project, "' ,@Year = ", year, "", sep=""), as.is=TRUE)
  odbcClose(channel)
  return(ds)}


