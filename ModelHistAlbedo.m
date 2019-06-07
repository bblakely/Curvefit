
%cd('/Users/bethanyblakely/Desktop/Analysis/AlbedoCurvefit')

DATAE=csvread('BinnedEvergreen_Albedo_fix.csv',1);
DATAD=csvread('BinnedDeciduous_Albedo_fix.csv',1);
DATAMx=csvread('BinnedMixed_Albedo_fix.csv',1);
DATAC=csvread('BinnedCrop_Albedo_fix.csv',1);
DATAMo=csvread('BinnedMosaic_Albedo_fix.csv',1);

%Read in Evergreeen LAI means manually
%NAME THE DATA VECTOR SOMETHING THATS NOT X
%Q=x;


LAIE=DATAE(:,2);
SNOWE=DATAE(:,3);
COUNTSE=DATAE(:,50);

LAID=DATAD(:,2);
SNOWD=DATAD(:,3);
COUNTSD=DATAD(:,50);

LAIMx=DATAMx(:,2);
SNOWMx=DATAMx(:,3);
COUNTSMx=DATAMx(:,50);

LAIC=DATAC(:,2);
SNOWC=DATAC(:,3);
COUNTSC=DATAC(:,50);


LAIMo=DATAMo(:,2);
SNOWMo=DATAMo(:,3);
COUNTSMo=DATAMo(:,50);

EvgreenfitsAlb=cell(46,1);
DecidfitsAlb=cell(46,1);
MixedfitsAlb=cell(46,1);
CropfitsAlb=cell(46,1);
MosaicfitsAlb=cell(46,1);

%Step one: Create fits

%Evergreen
SnowPlaceE=(1:46);
for i = 1:46
    
AlbindexE=DATAE(:,i+3);
Snowy=(AlbindexE(SNOWE>0));
SnowPlace=length(Snowy(~isnan(Snowy)));
SnowPlaceE(i)=SnowPlace;

if SnowPlace < 25 
     EvgreenfitsAlb{i}=EvgreenAlbFitGro(LAIE,AlbindexE);
elseif SnowPlace < 250     
     EvgreenfitsAlb{i}=EvgreenAlbFitWide(LAIE,SNOWE,AlbindexE, COUNTSE);
else
     EvgreenfitsAlb{i}=EvgreenAlbFit(LAIE,SNOWE,AlbindexE, COUNTSE);    

end
end


%Deciduous
SnowplaceD=(1:46);
for i = 1:46
    
AlbindexD=DATAD(:,i+3);
Snowy=(AlbindexD(SNOWD>0));
SnowPlace=length(Snowy(~isnan(Snowy)));
SnowPlaceD(i)=SnowPlace;

if SnowPlace < 25
      DecidfitsAlb{i}=DecidAlbFitGro(LAID,AlbindexD);
else
      DecidfitsAlb{i}=DecidAlbFit(LAID,SNOWD,AlbindexD, COUNTSD);    
end
end

%Mixed
SnowplaceMx=(1:46);
for i = 1:46
    
AlbindexMx=DATAMx(:,i+3);
Snowy=(AlbindexMx(SNOWMx>0));
SnowPlace=length(Snowy(~isnan(Snowy)));
SnowPlaceMx(i)=SnowPlace;

if SnowPlace < 100
      MixedfitsAlb{i}=MixedAlbFitGro(LAIMx,AlbindexMx);
elseif SnowPlace<400
      MixedfitsAlb{i}=MixedAlbFitWide(LAIMx,SNOWMx,AlbindexMx);
else 
        
     MixedfitsAlb{i}=MixedAlbFit(LAIMx,SNOWMx,AlbindexMx);    

end
end


%Crop
SnowplaceC=(1:46);
for i = 1:46
    
AlbindexC=DATAC(:,i+3);
Snowy=(AlbindexC(SNOWC>0));
SnowPlace=length(Snowy(~isnan(Snowy)));
SnowPlaceC(i)=SnowPlace;

if SnowPlace < 25
     CropfitsAlb{i}=CropAlbFitGro(LAIC,AlbindexC);
     
else 
     CropfitsAlb{i}=CropAlbFit(LAIC,SNOWC, AlbindexC);  

end
end


%Mosaic
SnowplaceMo=(1:46);
for i = 1:46
    
AlbindexMo=DATAMo(:,i+3);
Snowy=(AlbindexMo(SNOWMo>0));
SnowPlace=length(Snowy(~isnan(Snowy)));
SnowplaceMo(i)=SnowPlace;


if SnowPlace < 25
    MosaicfitsAlb{i}=MosaicAlbFitGro(LAIMo,AlbindexMo);
     
elseif SnowPlace <200
    MosaicfitsAlb{i}=MosaicAlbFitWide(LAIMo,SNOWMo,AlbindexMo); 
        
else       
    MosaicfitsAlb{i}=MosaicAlbFit(LAIMo,SNOWMo,AlbindexMo);  

end
end


%Step 2: assign values based on fits

%You must read in the paleoveg CSV for this section
nobs=length(B1);
FMW=zeros(nobs,46);

PaleoLat=Lat;
PaleoLon=Lon;

%Set these to wherever the snow and LAI files are; 
%Check overlap and dimensions for match prior to reding in;
PaleoLAIfile=csvread('F:\VegScale\CSVs_Alb\LAI_PNV.csv',1);
%('F:\VegScale\CSVs_Alb\LAI_v2_64k.csv',1);
PaleoSNOfile=csvread('F:\VegScale\CSVs_Alb\Snow_PNV.csv',1);
%('F:\VegScale\CSVs_Alb\Snow_64k.csv',1);

for x = 1:46
    DecidDate=DecidfitsAlb{x};
    EvgreenDate=EvgreenfitsAlb{x};
    MixedDate=MixedfitsAlb{x};
    CropDate=CropfitsAlb{x};
    MosaicDate=MosaicfitsAlb{x};
    
    PaleoLAI=PaleoLAIfile(:,x+6);
    PaleoSNO=PaleoSNOfile(:,x+6);
    
    x
for y = 1:nobs
    %Evergreen    
    if B1(y) == 1
        if SnowPlaceE(x)<25
            FMW(y,x)=EvgreenDate(PaleoLAI(y));
        else
            FMW(y,x)=EvgreenDate(PaleoLAI(y), PaleoSNO(y));
     
        end
    
    %Decid
    elseif B1(y)==4
        if SnowPlaceD(x)<25
            FMW(y,x)=DecidDate(PaleoLAI(y));
        else
            FMW(y,x)=DecidDate(PaleoLAI(y), PaleoSNO(y));
        end
    
    %Mixed
    elseif B1(y) == 5
        if SnowPlaceMx(x)<100
            FMW(y,x)=MixedDate(PaleoLAI(y));
        else
            FMW(y,x)=MixedDate(PaleoLAI(y), PaleoSNO(y));
            
        end
    
    
    %Crop
    elseif B1(y) == 12
        if SnowPlaceC(x)<25
            FMW(y,x)=CropDate(PaleoLAI(y));
        else
            FMW(y,x)=CropDate(PaleoLAI(y), PaleoSNO(y));
            
        end
    
   %Mosaic
     else if B1(y) == 14
        if SnowplaceMo(x)<25
            FMW(y,x)=MosaicDate(PaleoLAI(y));
        else
            FMW(y,x)=MosaicDate(PaleoLAI(y), PaleoSNO(y)); 
        end
        
    %otherwise
         else
         FMW(y,x)=9999;
         end
    end
    
end
            
        
end


FMW(isnan(FMW))=9999; %Necessary  because missing snow data results in some NaNs

%Set these each time; 

dim1=26
%18;%dim1 is columns or x dimension
dim2=20
%22;%dim2 is rows or x dimension

%These may or may not be obviously correct;
%UTM/GLSL will not have even degrees, WGS should.
%Check for smooth, slow acension all the way down columns for lat
%Or smooth, slow descension along rows for lon
LatCheck=reshape(PaleoLat,[dim1,dim2]);
LonCheck=reshape(PaleoLon, [dim1,dim2]);


%Step 3: print as gridded text files
%Don't forget to adjust filename if you want. 
%I don't, but then assign specific naming in the IDL reassembly script
for z = 1:46
FMWmat=reshape(FMW(:,z),[dim1,dim2]);
FMWmat(isnan(FMWmat))= 9999;
filename=strcat('C:\Users\Rocha Lab\Documents\MATLAB\Curvefit\Alb\Output\','Alb_PNV_',int2str(z), '.txt');
dlmwrite(filename,FMWmat)
end
    

