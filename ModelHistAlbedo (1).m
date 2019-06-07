
%cd('/Users/bethanyblakely/Desktop/Analysis/AlbedoCurvefit')

%You must replace 'NA's with NaN in excel before reading in new files

DATAE=csvread('BinnedEvergreen_AlbedoFlexRange.csv',1);
DATAD=csvread('BinnedDeciduous_AlbedoFlexRange.csv',1);
DATAMx=csvread('BinnedMixed_AlbedoFlexRange.csv',1);
DATAC=csvread('BinnedCrop_AlbedoFlexRange.csv',1);
DATAMo=csvread('BinnedMosaic_AlbedoFlexRange.csv',1);
DATAU=csvread('BinnedUrban_AlbedoFlexRange.csv',1);

%Read in Evergreeen LAI means manually
%NAME THE DATA VECTOR SOMETHING THATS NOT X
 LAI_EG=x2;

%Read in gridded datasets.
LAI_E=csvread('GridEvergreen_LAIFlexRange.csv',1,1);
SNOW_E=csvread('GridEvergreen_SnowFlexRange.csv',1,1);
COUNTS_E=csvread('GridEvergreen_CountsFlexRange.csv',1,1);
% LAIE=DATAE(:,2);
% SNOWE=DATAE(:,3);
% COUNTSE=DATAE(:,50);

LAI_D=csvread('GridDeciduous_LAIFlexRange.csv',1,1);
SNOW_D=csvread('GridDeciduous_SnowFlexRange.csv',1,1);
COUNTS_D=csvread('GridDeciduous_CountsFlexRange.csv',1,1);
% LAID=DATAD(:,2);
% SNOWD=DATAD(:,3);
% COUNTSD=DATAD(:,50);

LAI_Mx=csvread('GridMixed_LAIFlexRange.csv',1,1);
SNOW_Mx=csvread('GridMixed_SnowFlexRange.csv',1,1);
COUNTS_Mx=csvread('GridMixed_CountsFlexRange.csv',1,1);
% LAIMx=DATAMx(:,2);
% SNOWMx=DATAMx(:,3);
% COUNTSMx=DATAMx(:,50);

LAI_C=csvread('GridCrop_LAIFlexRange.csv',1,1);
SNOW_C=csvread('GridCrop_SnowFlexRange.csv',1,1);
COUNTS_C=csvread('GridCrop_CountsFlexRange.csv',1,1);
% LAIC=DATAC(:,2);
% SNOWC=DATAC(:,3);
% COUNTSC=DATAC(:,50);


LAI_Mo=csvread('GridMosaic_LAIFlexRange.csv',1,1);
SNOW_Mo=csvread('GridMosaic_SnowFlexRange.csv',1,1);
COUNTS_Mo=csvread('GridMosaic_CountsFlexRange.csv',1,1);
% LAIMo=DATAMo(:,2);
% SNOWMo=DATAMo(:,3);
% COUNTSMo=DATAMo(:,50);

LAI_U=csvread('GridUrban_LAIFlexRange.csv',1,1);
SNOW_U=csvread('GridUrban_SnowFlexRange.csv',1,1);
COUNTS_U=csvread('GridUrban_CountsFlexRange.csv',1,1);


EvgreenfitsAlb=cell(46,1);
DecidfitsAlb=cell(46,1);
MixedfitsAlb=cell(46,1);
CropfitsAlb=cell(46,1);
MosaicfitsAlb=cell(46,1);
UrbanfitsAlb=cell(46,1);

%For some reason you must declare r2st global before running this loop
global r2Alb

%Curvefit loops. Use binned data to create 3d fits of snow and LAI to
%albedo by vegytype. Takes a long time (several minutes)

%Evergreen
SnowplaceE=(1:46);
FitStatsE=cell(46,1);

%Specific fit metrics
sseE=(1:46);
r2E=(1:46);
rmseE=(1:46);

for i = 1:46

SNOWE=SNOW_E(:,i); 
LAIE=LAI_E(:,i);

AlbindexE=DATAE(:,i+1);
Snowy=(AlbindexE(SNOWE>0));
SnowPlace=length(Snowy(~isnan(Snowy)));
SnowPlaceE(i)=SnowPlace;

COUNTSE=COUNTS_E(:,i);

if SnowPlace < 25 
     EvgreenfitsAlb{i}=EvgreenAlbFitGro(LAIE,AlbindexE);
elseif SnowPlace < 250     
     EvgreenfitsAlb{i}=EvgreenAlbFitWide(LAIE,SNOWE,AlbindexE, COUNTSE);
else
     EvgreenfitsAlb{i}=EvgreenAlbFit(LAIE,SNOWE,AlbindexE, COUNTSE);    

end

%Save fit stats
FitStatsE{i}=r2Alb;

sseE(i)=r2Alb.sse;
r2E(i)=r2Alb.rsquare;
rmseE(i)=r2Alb.rmse;

end



%Deciduous
SnowplaceD=(1:46);
FitStatsD=cell(46,1);

sseD=(1:46);
r2D=(1:46);
rmseD=(1:46);
for i = 1:46

SNOWD=SNOW_D(:,i); 
LAID=LAI_D(:,i);

AlbindexD=DATAD(:,i+1);
Snowy=(AlbindexD(SNOWD>0));
SnowPlace=length(Snowy(~isnan(Snowy)));
SnowPlaceD(i)=SnowPlace;

COUNTSD=COUNTS_D(:,i);

if SnowPlace < 25
      DecidfitsAlb{i}=DecidAlbFitGro(LAID,AlbindexD);
 elseif SnowPlace < 250     
      DecidfitsAlb{i}=DecidAlbFitWide(LAID,SNOWD,AlbindexD, COUNTSD);
else
      DecidfitsAlb{i}=DecidAlbFit(LAID,SNOWD,AlbindexD,COUNTSD);    
end

%Save fit stats
FitStatsD{i}=r2Alb;

sseD(i)=r2Alb.sse;
r2D(i)=r2Alb.rsquare;
rmseD(i)=r2Alb.rmse;

end

%Mixed
SnowplaceMx=(1:46);
FitStatsMx=cell(46,1);

sseMx=(1:46);
r2Mx=(1:46);
rmseMx=(1:46);
for i = 1:46
    
SNOWMx=SNOW_Mx(:,i); 
LAIMx=LAI_Mx(:,i);

AlbindexMx=DATAMx(:,i+1);
Snowy=(AlbindexMx(SNOWMx>0));
SnowPlace=length(Snowy(~isnan(Snowy)));
SnowPlaceMx(i)=SnowPlace;

COUNTSMx=COUNTS_Mx(:,i);

if SnowPlace < 25
      MixedfitsAlb{i}=MixedAlbFitGro(LAIMx,AlbindexMx);
elseif SnowPlace<250
      MixedfitsAlb{i}=MixedAlbFitWide(LAIMx,SNOWMx,AlbindexMx,COUNTSMx);
else 
        
     MixedfitsAlb{i}=MixedAlbFit(LAIMx,SNOWMx,AlbindexMx,COUNTSMx);    

end
FitStatsMx{i}=r2Alb;

sseMx(i)=r2Alb.sse;
r2Mx(i)=r2Alb.rsquare;
rmseMx(i)=r2Alb.rmse;
end


%Crop
SnowplaceC=(1:46);
FitStatsC=cell(46,1);

sseC=(1:46);
r2C=(1:46);
rmseC=(1:46);
for i = 1:46
    
SNOWC=SNOW_C(:,i); 
LAIC=LAI_C(:,i);

AlbindexC=DATAC(:,i+1);
Snowy=(AlbindexC(SNOWC>0));
SnowPlace=length(Snowy(~isnan(Snowy)));
SnowPlaceC(i)=SnowPlace;

COUNTSC=COUNTS_C(:,i);

if SnowPlace < 25
     CropfitsAlb{i}=CropAlbFitGro(LAIC,AlbindexC);
     
elseif SnowPlace<250
    CropfitsAlb{i}=CropAlbFitWide(LAIC,SNOWC,AlbindexC,COUNTSC);
     
else 
     CropfitsAlb{i}=CropAlbFit(LAIC,SNOWC, AlbindexC, COUNTSC);  

end
FitStatsC{i}=r2Alb;

sseC(i)=r2Alb.sse;
r2C(i)=r2Alb.rsquare;
rmseC(i)=r2Alb.rmse;
end


%Mosaic
SnowplaceMo=(1:46);
FitStatsMo=cell(46,1);

sseMo=(1:46);
r2Mo=(1:46);
rmseMo=(1:46);

for i = 1:46
    
SNOWMo=SNOW_Mo(:,i); 
LAIMo=LAI_Mo(:,i);
    
AlbindexMo=DATAMo(:,i+1);
Snowy=(AlbindexMo(SNOWMo>0));
SnowPlace=length(Snowy(~isnan(Snowy)));
SnowplaceMo(i)=SnowPlace;

COUNTSMo=COUNTS_Mo(:,i);

if SnowPlace < 25
    MosaicfitsAlb{i}=MosaicAlbFitGro(LAIMo,AlbindexMo);
     
elseif SnowPlace <250
    MosaicfitsAlb{i}=MosaicAlbFitWide(LAIMo,SNOWMo,AlbindexMo, COUNTSMo); 
        
else       
    MosaicfitsAlb{i}=MosaicAlbFit(LAIMo,SNOWMo,AlbindexMo, COUNTSMo);  

end
FitStatsMo{i}=r2Alb;

sseMo(i)=r2Alb.sse;
r2Mo(i)=r2Alb.rsquare;
rmseMo(i)=r2Alb.rmse;
end


%Urban
SnowplaceU=(1:46);
FitStatsU=cell(46,1);

sseU=(1:46);
r2U=(1:46);
rmseU=(1:46);
for i = 1:46
    
SNOWU=SNOW_U(:,i); 
LAIU=LAI_U(:,i);
    
AlbindexU=DATAU(:,i+1);
Snowy=(AlbindexU(SNOWU>0));
SnowPlace=length(Snowy(~isnan(Snowy)));
SnowplaceU(i)=SnowPlace;

COUNTSU=COUNTS_U(:,i);

if SnowPlace < 25
    UrbanfitsAlb{i}=UrbanAlbFitGro(LAIU,AlbindexU);
     
elseif SnowPlace <250
    UrbanfitsAlb{i}=UrbanAlbFitWide(LAIU,SNOWU,AlbindexU, COUNTSU); 
        
else       
    UrbanfitsAlb{i}=UrbanAlbFit(LAIU,SNOWU,AlbindexU, COUNTSU);  

end
FitStatsU{i}=r2Alb;

sseU(i)=r2Alb.sse;
r2U(i)=r2Alb.rsquare;
rmseU(i)=r2Alb.rmse;
end


%You must read in the paleoveg CSV for this section
nobs=length(B1);
FMW=zeros(nobs,46);

PaleoLat=Lat;
PaleoLon=Lon;


ModernLAIfile=csvread('ModernLAI_clean.csv',1);
PaleoLAIfile=csvread('Paleo_LAI_UTM_v2.csv',1); %clean1

PaleoSNOfile=csvread('SNOW_clean_1.csv',1);
ModernSnowFile=PaleoSNOfile;

%Assignment loop. Inputs LAI, snow, and vegtype from files above to
%curvefits created in the prevous loop.

for x = 1:46
    DecidDate=DecidfitsAlb{x};
    EvgreenDate=EvgreenfitsAlb{x};
    MixedDate=MixedfitsAlb{x};
    CropDate=CropfitsAlb{x};
    MosaicDate=MosaicfitsAlb{x};
    UrbanDate=UrbanfitsAlb{x};
    
    PaleoLAI= PaleoLAIfile(:,x+2);
    %ModernLAIfile(:,x+2);
    %PaleoLAIfile(:,x+2);
    PaleoSNO= PaleoSNOfile(:,x+2);
    %PaleoSNOfile(:,x+2);
   
    x
for y = 1:nobs
    %Evergreen    
    if B1(y) == 1
        if SnowPlaceE(x)<25
            FMW(y,x)=EvgreenDate(LAI_EG(x));
        elseif SnowPlaceE(x)<250
            %PaleoLAI(y)
            FMW(y,x)=EvgreenDate (LAI_EG(x), PaleoSNO(y));
        else
            FMW(y,x)=EvgreenDate(LAI_EG(x), PaleoSNO(y));
     
        end
    
    %Decid
    elseif B1(y)==4
        if SnowPlaceD(x)<25
            FMW(y,x)=DecidDate(PaleoLAI(y));
         elseif SnowPlaceE(x)<250
             FMW(y,x)=DecidDate (PaleoLAI(y), PaleoSNO(y));
        else
            FMW(y,x)=DecidDate(PaleoLAI(y), PaleoSNO(y));
        end
    
    %Mixed
    elseif B1(y) == 5
        if SnowPlaceMx(x)<25
            FMW(y,x)=MixedDate(PaleoLAI(y));
        elseif SnowPlaceE(x)<250
            FMW(y,x)=MixedDate(PaleoLAI(y), PaleoSNO(y));
        else
            FMW(y,x)=MixedDate(PaleoLAI(y), PaleoSNO(y));
            
        end
    
    
    %Crop
    elseif B1(y) == 12
        if SnowPlaceC(x)<25
            FMW(y,x)=CropDate(PaleoLAI(y));
        elseif SnowPlaceE(x)<250
            FMW(y,x)=CropDate(PaleoLAI(y), PaleoSNO(y));
        else
            FMW(y,x)=CropDate(PaleoLAI(y), PaleoSNO(y));
            
        end
    
   %Mosaic
    elseif B1(y) == 14
        if SnowplaceMo(x)<25
            FMW(y,x)=MosaicDate(PaleoLAI(y));
        elseif SnowPlaceE(x)<250
            FMW(y,x)=MosaicDate(PaleoLAI(y), PaleoSNO(y));
        else
            FMW(y,x)=MosaicDate(PaleoLAI(y), PaleoSNO(y)); 
        end
        
    %Urban
    elseif B1(y) == 13
        if SnowplaceU(x)<25
            FMW(y,x)=UrbanDate(PaleoLAI(y));
        elseif SnowPlaceE(x)<250
            FMW(y,x)=UrbanDate(PaleoLAI(y), PaleoSNO(y));
        else
            FMW(y,x)=UrbanDate(PaleoLAI(y), PaleoSNO(y)); 
        end
        
    %otherwise
         else
         FMW(y,x)=9999;
         end
    end
    
end
            
        
%end



LatCheck=reshape(PaleoLat,[141,174]);
LonCheck=reshape(PaleoLon, [141,174]);

%Writing loop. Must make folder to receive files
for z = 1:46
FMWmat=reshape(FMW(:,z),[141,174]);
FMWmat(isnan(FMWmat))= 9999;
filename=strcat('Output/','Albdo_paleo_v2_',int2str(z), '.txt')
%'/Users/bethanyblakely/Desktop/Analysis/AlbedoPaleo_v2/','Albdo_paleo_v2_',int2str(z), '.txt')
dlmwrite(filename,FMWmat)
end

%Knitting done in IDL: read.ascii.pro


%Summary of metrics
min(rmseC)
find(rmseC==min(rmseC))
max(rmseC)
find(rmseC==max(rmseC))

mean(rmseC)
std(rmseC)

min(rmseU)
find(rmseU==min(rmseU))
max(rmseU)
find(rmseU==max(rmseU))

mean(rmseU)
std(rmseU)

min(rmseMo)
find(rmseMo==min(rmseMo))
max(rmseMo)
find(rmseMo==max(rmseMo))

mean(rmseMo)
std(rmseMo)

min(rmseMx)
find(rmseMx==min(rmseMx))
max(rmseMx)
find(rmseMx==max(rmseMx))

mean(rmseMx)
std(rmseMx)

