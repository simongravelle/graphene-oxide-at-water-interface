clear all
close all

% RT=0.593; % kcal/mol
T=300;
kB=3.297623483e-24; % cal/K
kBT=kB*T/1000*6.022e23; % kCal/mol

N0=45;
allOC=[0,5,17,34,55];
allOH=[10,20,30,40,50,60,70];
allM=['a','b','c','d','e','f','g','h'];
allL=[1.6,2,2.4,2.8];

minZ0=-28;
maxZ0=28;

cpt=0;
for ii=minZ0:2:maxZ0
	cpt=cpt+1;
	if ii<0
		DataN{cpt}=['Datam',num2str(-ii)];
		Z0{cpt}=ii;
	else
		DataN{cpt}=['Data',num2str(ii)];
		Z0{cpt}=ii;
	end
end

cutoff=500;

for ll=1:length(allL)
	L=allL(ll);
	for oc=1:length(allOC)
		OC=allOC(oc);
		for oh=1:length(allOH)
			OH=allOH(oh);
			for m=1:length(allM)
				M=allM(m);
				%folder=['L',num2str(L),'_OC',num2str(OC),'_OH',num2str(OH),'dis_',M];
				folder=['L',num2str(L),'_OC',num2str(OC),'_OH',num2str(OH),M];
				X=['../',folder,'/Datam',num2str(28),'/Result.dat'];
				if exist(X,'file')==2
					file=fopen('metadata.dat','wt');
					clear allpro;
					clear avdenpro;
					for aa=1:length(DataN)
						Subfolder=DataN{aa};
						% isolate the correct portion of the data
						X=['../',folder,'/',Subfolder,'/Result.dat'];
						dataOriginal=load(X);
						dataTreated=dataOriginal([cutoff:end],[1:2]);
						X=['../',folder,'/',Subfolder,'/resultTreated.dat'];
						save([X],'dataTreated','-ascii');
						% update the megadata file
						z0=Z0{aa};
						X=['../',folder,'/',Subfolder,'/resultTreated.dat ',num2str(z0),' 5'];
						fprintf(file,X);
						fprintf(file,'\n');
						% density profile
						X=['../',folder,'/',Subfolder,'/waterprofile.dat'];
						denpro=dlmread(X,' ',4,1);
						if length(denpro)>0
							denpro(:,1)=[];
							denpro(:,1)=[];
							denpro(:,2)=[];
							xcoor=denpro(:,1);
							allpro(:,aa)=denpro(:,2);
						end
					end
					X=['../',folder,'/waterprofile.dat'];
					avdenpro(:,1)=xcoor;
					avdenpro(:,2)=mean(allpro')';
					save([X],'avdenpro','-ascii');
					fclose(file);
					% run the wham algorithm
					X=['./wham ',num2str(minZ0),' ',num2str(maxZ0),' ',num2str(N0),' 1e-8 300 0 metadata.dat PMF',num2str(N0),'.dat 5 455878'];
					system(X);
					% copy the PMF
					Y=['PMF',num2str(N0),'.dat'];
					X=['../',folder,'/'];
					copyfile(Y,X);
				end
			end
		end
	end
end
		
for ll=1:length(allL)
	L=allL(ll);
	for oc=1:length(allOC)
		OC=allOC(oc);
		for oh=1:length(allOH)
			OH=allOH(oh);
			clear AvPMF
			clear tempPMF
			clear Avzsurf;
			clear DEAds;
			n=0;
			for m=1:length(allM)
				M=allM(m);
				%folder=['L',num2str(L),'_OC',num2str(OC),'_OH',num2str(OH),'dis_',M];
				folder=['L',num2str(L),'_OC',num2str(OC),'_OH',num2str(OH),M];
				X=['../',folder,'/PMF',num2str(N0),'.dat'];
				Y=['../',folder,'/waterprofile.dat'];
				
				if exist(X,'file')==2 && exist(Y,'file')==2
					
					% detect the interface zsurf
					denpro=load(Y);
					denprocenter=denpro;
					denprocenter(denprocenter(:,1)<-15,:)=[];
					denprocenter(denprocenter(:,1)>15,:)=[];				
					halfden=denprocenter(2,2)/2;
					intpro(:,1)=[min(denprocenter(:,1)):0.1:max(denprocenter(:,1))]';
					for iji=1:length(intpro(:,1))
						intpro(iji,2)=interp1(denprocenter(:,1),denprocenter(:,2),intpro(iji,1));
					end
					dist=sqrt((intpro(:,2)-halfden).^2);
					[minv,idx]=min(dist);
					zsurf=intpro(idx,1);					
				
					PMF=load(X);
					n=n+1;
					AvPMF(:,1)=PMF(:,1);
					tempPMF(:,n)=PMF(:,2);
					Avzsurf(n)=zsurf;
					AvDens(:,1)=denpro(:,1);
					tempDens(:,n)=denpro(:,2);
				end
			end
			
			if n>0			
				AvPMF(:,1)=AvPMF(:,1)-mean(Avzsurf);
				AvPMF(:,2)=mean(tempPMF')'/kBT; % kbT units
				AvPMF(:,3)=std(tempPMF')'/kBT; % kbT units
				
				%X=['../Data_dis/PMF_L',num2str(L),'_OC',num2str(OC),'_OH',num2str(OH),'.dat'];
				X=['../Data/PMF_L',num2str(L),'_OC',num2str(OC),'_OH',num2str(OH),'.dat'];
				save([X],'AvPMF','-ascii');
			
				AvDens(:,1)=AvDens(:,1)-mean(Avzsurf);
				AvDens(:,2)=mean(tempDens')';
				AvDens(:,3)=std(tempDens')';
				
				%X=['../Data_dis/DensityPro_L',num2str(L),'_OC',num2str(OC),'_OH',num2str(OH),'.dat'];
				X=['../Data/DensityPro_L',num2str(L),'_OC',num2str(OC),'_OH',num2str(OH),'.dat'];
				save([X],'AvDens','-ascii');
			
			end
		end
	end
end
				
				
				
				
				
				
				
				
