# NHANES_wtofm

The following procedures and supporting files are needed to successfully run a survey-weighted overfitted latent class model using 24-hour dietary intake data from the National Health and Nutrition Examination Surveys, pooling four survey cycles: 2011-2012, 2013-2014, 2015-2016, 2017-2018.

## STEP 1 - Download FPED data
FPED data can be downloaded directly from the [USDA/ARS](https://www.ars.usda.gov/northeast-area/beltsville-md-bhnrc/beltsville-human-nutrition-research-center/food-surveys-research-group/docs/fped-databases/ "USDA/ARS title") website.
You will need all the food pattern equivalents for foods in the WWEIA, NHANES 20XX-YY. Example attached pools 2011-2018 and can be found in [nhanes-data](www.github.com/bjks10/NHANES_wtofm/nhanes-data "nhanes-data title") subfolder. 

## STEP 2 - Make Food Groups from 24HR Diet Recall data
Run MakeFoodGroups_nhanes.sas  - run for each survey cycle separately. Make sure the order of variables match with the cycle run. Example runs 2011-12 (G), 2013-14 (H), 2015-16 (I), 2017-18 (J)can be found in [nhanes-data](www.github.com/bjks10/NHANES_wtofm/nhanes-data "nhanes-data title") subfolder.
Note: The ordering of food group variables changed in 2017-2018 (J) cycle
	
  *Input files*:	FPED_DR1TOT_XXYY.sas7bdat
					        FPED_DR2TOT_XXYY.sas7bdat
					        FPED_DR1IFF_XXYY.sas7bdat
					        FPED_DR2IFF_XXYY.sas7bdat

	*Output files*:	FPED_DRAVGXXYY.sas7bdat

## STEP 3 - Create Tertiles of Consumption
Run fped_CreateTerts.sas - Combine survey cycles to a pooled dataset, drop "total" food labels, add diet weights. 
XXYY denote a single survey cycle (e.g. 1112 for 2011-2012 cycle). 
XXZZ is the full range of combined surveys (e.g. 1118 for 2011-2018 surveys pooled). 
[A] denotes the alpha code attached to each survey cycle: G=2011-2012, H=2013-2014, I=2015-2016, J=2017-2018.

	*Input files*: 	DEMO_[A].sas7bdat
				          FPED_DRAVGXXYY.sas7bdat
				          DR1TOT_[A].sas7bdat

	*Output files*:	tert_nhanesXXZZ.sas7bdat

## STEP 4 - Calculate HEI2015 scores for all participants
Run GenerateHEI.sas - calculate HEI2015 score for each survey cycle [X].
	
  *Input files*:	DEMO_X.sas7bdat
				          FPED_DR1TOT_X.sas7bdat
				          FPED_DR2TOT_X.sas7bdat
				          hei2015_score_macro.sas

	*Output files*:	HEIXXYY_avg.sas7bdat

## STEP 5 - Merge HEIscores to FPED tertile data
Run MergeFPED_HEI.sas - Merge FPED-tertiles and HEI2015 score and demographic data to single dataset and export to CSV file
	
  *Input files*:	HEIXXYY_avg.sas7bdat
				          tert_nhanesXXZZ.sas7bdat

	*Output files*: hei_tert1118.sas7bdat
				          adultHEI_fped1118terts.csv

## STEP 6 - Run Survey-weighted Overfitted Latent Class model
Import csv file into MATLAB and run wtOFM_nhanes.m under K=50 clusters.

	*Input files*:	adult HEI_fped1118terts.csv
	*Output files*:	nhanesadult_50_wtOFMresults.xlsx
				          wtOFM_50_NHANESadultresults.mat
		              diet[#]dist.png (for each nonempty cluster profile derived)
				          noconsum_adultpat.png
				          highconsum_adultpat.png
				          theta0_adultpattern.pdf
