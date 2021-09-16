%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Robust Profile Clustering 
%% Programmer: BJKS     
%% Data: NHANES - all adults       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc;
% load('sim_wRPCdataRand')
adultHEIFPED=readtable('adultHEI_fped1118terts1Sep2021.csv');
%load('nhanes_wtFPED')

low_income=adultHEIFPED.povbin;
SEQN=adultHEIFPED.SEQN;
dietwt=adultHEIFPED.dietwt8yr;
food=table2array(adultHEIFPED(:,21:49));
gender=adultHEIFPED.RIAGENDR;
age=adultHEIFPED.RIDAGEYR;
intwt=adultHEIFPED.wtint8yr;
all_indata=[SEQN dietwt low_income gender age food];



%normalization constant
food=all_indata(:,6:end);
[n,p]=size(food);


wtsum=sum(dietwt);
N=235919796;
c=wtsum/N;

wt_c=dietwt/c;
    k_max=50;
    d_max=max(food(:));
    d=max(food);

wt_theta0=repmat(wt_c,[1,p]);
    %vectorization of data
     idz = repmat(1:p,n,1); idz = idz(:);
     y_d = food(:); lin_idx = sub2ind([p,d_max],idz,y_d);

    %% SET UP PRIORS %%

        %pi_h for all classes
        sp_k=50;
a_pi=ones(1,k_max)/sp_k;
pi_h=drchrnd(a_pi,1);

        %phi - cluster index

        rr = unifrnd(0,1,[n,1]);
       pisums=[0 cumsum(pi_h)];
       C_i=zeros(n,1);
x_ci=zeros(n,k_max);
        for l = 1:k_max
            ind = rr>pisums(l) & rr<=pisums(l+1);
            C_i(ind==1) = l;
            x_ci(:,l)=ind;
        end 
n_C_i=sum(x_ci);


          %global theta0/1
     eta=ones(1,d_max);
    theta0=zeros(p,k_max,d_max);

    for k=1:k_max
        for j=1:p
            dj=d(j);
            theta0(j,k,1:dj)=drchrnd(eta(1:dj),1);
         end
    end

 
    %% ------------ %%
    %% data storage %%
    %% ------------ %%
    nrun=25000; burn=15000;  thin=5;
    pi_out=zeros(nrun/thin,k_max);
   
    theta0_out=zeros(nrun/thin,p,k_max,d_max);
    ci_out=zeros(nrun/thin,n);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% POSTERIOR COMPUTATION %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    As=zeros(n,p);
    p_ij=zeros(n,p);
tic
for iter=1:nrun

   
        %% -- update pi_h -- %%
         for h=1:k_max
%              n_C_i(h)=sum(C_i==h);
            n_C_i(h)=sum(wt_c(C_i==h));
         end

        a_pih=a_pi+n_C_i;
        pi_h=drchrnd(a_pih,1);

      %% -- phi ~multinomial(pi_h) -- %%

  Cp_k=zeros(n,k_max);
for k=1:k_max
    t0h=reshape(theta0(:,k,:),p,d_max);
    tmpmat0=reshape(t0h(lin_idx),[n,p]);
    Cp_k(:,k)=pi_h(k)*prod(tmpmat0,2);
end
probCi = bsxfun(@times,Cp_k,1./(sum(Cp_k,2)));
    x_ci=mnrnd(1,probCi); [r, c]=find(x_ci); x_gc=[r c];
    x_gc=sortrows(x_gc,1); C_i=x_gc(:,2);

%store global cluster index for postprocess/relabelling
        % - update theta - %
    dmat0=zeros(p,d_max);
    for k=1:k_max
        C_is=repmat(C_i,[1,p]);
%          ph0 = (C_is==k); %subj's in global cluster h
           ph0=(C_is==k).*wt_theta0; %global cluster with weights
            for c = 1:d_max
                 dmat0(:,c) = sum((food==c).*ph0)';
            end
            for j=1:p
                dj=d(j);
                a_tn0=eta(1:dj)+dmat0(j,1:dj);
                theta0(j,k,1:dj) = drchrnd(a_tn0,1);
            end
    end

      
      if mod(iter,thin)==0
        pi_out(iter/thin,:)=pi_h;
        ci_out(iter/thin,:)=C_i;
        theta0_out(iter/thin,:,1:size(theta0,2),:)=theta0;
       end


%% RELABELLING STEP TO ENCOURAGE MIXING %%
    if mod(iter,10)==0
        new_order=randperm(k_max);
        newC_i=C_i;
        
        for k=1:k_max
            newC_i(C_i==k)=new_order(k);
        end
        
        C_i=newC_i;
        theta0=theta0(:,new_order,:);
    end




end
eltime=toc;
    pi_burn=pi_out((burn/thin)+1:end,:);
    theta0_burn=theta0_out((burn/thin)+1:end,:,:,:);
    ci_burn=ci_out((burn/thin)+1:end,:);

  
        figure; %check mixing of pi parameter
        plot(pi_burn)
%         saveas(gcf,'wtnhaneslow_pis.png')



 save('wtOFM_nhanesadult_MCMCout','pi_burn','ci_burn','theta0_burn','eltime','-v7.3');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% POST-PROCESS: PASAPILIOPOULIS SWITCH %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m=size(pi_burn,1); m_perm=size(theta0_burn,1); m_thin=m/m_perm;
k_med=median(sum(pi_burn>0.05,2));
pd=pdist(transpose(ci_burn),'hamming'); %prcnt differ
cdiff=squareform(pd); %Dij
Zci=linkage(cdiff,'complete');

figure; %save dendrogram of hierarchical clustering
dendrogram(Zci);
% saveas(gcf,'wtnhaneslow_dendrogram.png')

k0=k_med;
clust0 = cluster(Zci,'maxclust',k_med); %choose k0=5;


%Ordered MCMC set
ci_relabel=zeros(m,k0); %ci_relabelin=zeros(m,k_in); 
for l=1:k0
    ci_relabel(:,l) = mode(ci_burn(:,clust0==l),2);
end

pi_order=zeros(m,k0);
theta0_order=zeros(m,p,k0,d_max);
for iter=1:m
    iter_order=ci_relabel(iter,:);
    pi_h1=pi_burn(iter,iter_order);
    pi_order(iter,:)=pi_h1/sum(pi_h1);
    theta0_order(iter,:,:,:)=theta0_burn(iter,:,iter_order,:);
    
end

theta0_med=reshape(median(theta0_order),[p,k0,d_max]);
pi_med=median(pi_order);
% pi_med=pi_med/sum(pi_med);

[p,k0,d]=size(theta0_med);


%%Modal pattern of Global Clusters %%
[val0,ind0]=max(theta0_med,[],3);
t_ind0=transpose(ind0);
[t0,ia,ic]=unique(t_ind0,'rows');
k0=length(ia);

pi_med=pi_med(ia)/sum(pi_med(ia));
theta0_med=theta0_med(:,ia,:);
theta0_medi=theta0_med./sum(theta0_med,3);

theta0_probs=cell(d,1);
for dk=1:d
    theta0_probs{dk}=reshape(theta0_medi(:,:,dk),[p,k0]);
end
profile_patts = cell(k0,1);
for k=1:k0
   profile_patts{k}=reshape(theta0_medi(:,k,:),[p,d]);
end

[val_theta,pat_theta]=max(theta0_medi,[],3);

flabels={'Citrus, Melon Berries','Other fruit','Fruit juice','Dk Green Veg',...
    'Tomatoes','Other Red/Org veg','Potatoes','Other starchy veg','Other veg',...
    'Legumes (veg)','Whole grains','Refined grains','Meat(ns)','Cured meats',...
    'Organ meat','Poultry','Seafood (highn3)','Seafood (lown3)','Eggs','Soybean',...
    'Nuts/seeds','Legumes (protein)','Milk','Yogurt','Cheese','Oils','Solid fat',...
    'Added sugar','Alcohol'};
for k=1:k0
    figure; barh(profile_patts{k},'stacked')
    xlim(0:1)
    yticks(1:29)
    yticklabels(flabels)
    title(strcat('Dietary Profile  ',num2str(k)))
    xlabel('Posterior Probability')
    legend({'None','Low','Med','High'},'Location','southwestoutside')
    saveas(gcf,strcat('diet',num2str(k),'dist.png'))
end

  Cp_med=zeros(n,k0);
for k=1:k0
    t0h=reshape(theta0_med(:,k,:),p,d_max);
    tmpmat0=reshape(t0h(lin_idx),[n,p]);
    Cp_med(:,k)=pi_med(k)*prod(tmpmat0,2);
end
pCi = bsxfun(@times,Cp_med,1./(sum(Cp_med,2)));
[z_val,z_max]=max(pCi,[],2);
    x_ci=mnrnd(1,pCi); [r, c]=find(x_ci); x_gc=[r c];
    x_gc=sortrows(x_gc,1); z_med=x_gc(:,2);

nhanesfped_wofm = table(SEQN,dietwt,age,gender,z_med);
m1z=z_med;



filename = strcat('nhanesadult_',num2str(sp_k),'_wtOFMcresults.xlsx');

writetable(nhanesfped_wofm,filename)

save(strcat('wtOFM_', num2str(k_max),'_NHANESadultResults'),'dietwt','nhanesid','age','gender','m1z','pi_med','theta0_med', 'val_theta','pat_theta');
%% PLOTS %%

% Posterior probability plot of no consumption
figure; plot(theta0_probs{1},'Linewidth',1)
xticks(1:29)
xticklabels(flabels)
ylabel('Posterior Probability of No consumption')
saveas(gcf,'noconsum_adultpat.png')

% Posterior probability plot of high consumption
figure; plot(theta0_probs{4},'Linewidth',1)
xticks(1:29)
xticklabels(flabels)
ylabel('Posterior Probability of High consumption')
saveas(gcf,'highconsum_adultpat.png')


% Heatmap of theta0 - global pattern mode

figure; 
    h=heatmap(pat_theta)
    h.YDisplayLabels = flabels;
    h.XLabel = "Dietary Profile";
    h.Colormap = parula
saveas(gcf,'theta0_adultpattern.pdf')
