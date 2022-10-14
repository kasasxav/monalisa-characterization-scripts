datadir = ['\\storage.ad.scilifelab.se\monalisa\Big_MoNaLISA\2021project\Project\resolution\3-images-data\2020-08-05\original-data\cell10'];
outdir = [datadir '\fit'];
mkdir(outdir)


for z = 1
    
    fn = ['\line-profiles-2-r.csv']; %lista(1).name(1:end-4);%'02_pwr50';
%     fnC = ['cc_lines_slice' num2str(z) '.txt'];
%     fnD = ['dec_lines_slice' num2str(z) '.txt'];
    ll = importdata([datadir '\' fn]);
%     llC = importdata([datadir '\' fnC]);
%     llD = importdata([datadir '\' fnD]);
    pxSize = 25;
    confidence = 2.5;
    z = 1;
    k = 1;
    for i = 2:size(ll.data,2)
        %%
        profile = ll.data(1:end,i); 
        %profile = profile(profile>0);
        profile = (profile-min(profile))./max((profile-min(profile)));
%         profileC = llC.data(1:end,i); 
%         profileC = profileC(profile>0);  profileC = (profileC-min(profileC))./max((profileC-min(profileC)));
%         profileD = llD.data(1:end,i); 
%         profileD = profileD(profile>0);  profileD = (profileD-min(profileD))./max((profileD-min(profileD)));
%         if isempty(profile) || length(profile) < length(ll.data(1:end,i))*0.5
%             disp(i);
%             param(k,:,z) = zeros(1,10);
%             k=k+1;
%             continue
%         else
            nm = [1:length(profile)]*pxSize; nm=nm';
            %nmC = [1:length(profileC)]*pxSize; nmC=nmC';
            [maxima, maxId] = findpeaks(profile,'MinPeakHeight',max(profile)*0.025,'MinPeakDistance',(length(nm)-5));
            

            bkg = mean([profile(1:3); profile(end-3:end)]);
%             try
                % Single Lorentzian fit  p = a, w, x0, y0
                            ft1L = fittype( 'y0+(2*a/pi)*(w./(4*(x-x0).^2 + w.^2))', 'independent', 'x', 'dependent', 'y' );
                            opts1L = fitoptions( 'Method', 'NonlinearLeastSquares' );
                            opts.Algorithm = 'Levenberg-Marquardt';%'Trust-Region';%
                            opts1L.Display = 'Off';
                            opts1L.Lower = [0 48 0 0];
                            opts1L.Upper = [Inf 3000 length(profile)*pxSize Inf];
                            opts1L.MaxIter = 1000;
                            opts1L.Robust = 'Bisquare';
                            opts1L.StartPoint = [max(profile)*100 pxSize*2 maxId*pxSize bkg];
                            % Fit model to data.
                            [fitresult, gof] = fit( nm, profile, ft1L, opts1L );
                            profile_1L = fitresult(nm);
                            coeffvals_1L = coeffvalues(fitresult);
                
                
                % 2 gaussian fit  p = (a,f1,w1,w2,x0,y0,x)  =====
                %                      1   2  3  4  5   6
%                 ft2G = fittype( 'y0 + a*((1-f1)*exp(-((x-x0).^2)./(2*(w1.^2))) + (f1)*exp(-((x-x0).^2)./(2*(w2.^2))) )', 'independent', 'x', 'dependent', 'y' );
%                 opts2G = fitoptions( 'Method', 'NonlinearLeastSquares' );
%                 opts2G.Algorithm = 'Trust-Region'; 
%                 opts2G.Weights = (profile./max(profile)).^2;
%                 opts2G.Display = 'Off';
%                 opts2G.MaxIter = 1000;
%                 opts2G.Robust = 'LAR';
%                 opts2G.Lower = [max(profile)/2 0 pxSize/2.35 350/2.35  0 0];
%                 opts2G.Upper = [max(profile)*2 0.2 400/2.35 400/2.35 length(profile)*pxSize 0.15]; %min([profile(1:3) ; profile(end-3:end)])
%                 opts2G.StartPoint = [max(profile) 0.1 50/2.35 350/2.35 maxId*pxSize min(profile)]; 
%                 % Fit model to data.
%                 [fitresult2G, gofLC] = fit( nm, profile, ft2G, opts2G );
%                 profile_2G = fitresult2G(nm);
%                 
%                 fwhm(1:2) = sqrt(log(4)) * 2 *[fitresult2G.w1 fitresult2G.w2] ;
%                 fitG1 = fitresult2G.y0 + fitresult2G.a*((1-fitresult2G.f1)*exp(-((nm-fitresult2G.x0).^2)./(2*(fitresult2G.w1.^2))));
%                 fitG2 = fitresult2G.y0 + fitresult2G.a*(fitresult2G.f1*exp(-((nm-fitresult2G.x0).^2)./(2*(fitresult2G.w2.^2))));
%                 
                % Estimation of the bkg
%                 cuttedInt = round(fitresult2G.w2*confidence/(2*pxSize));
%                 bkg = [profile(1:maxId-cuttedInt); profile(maxId+cuttedInt:end)];
%                 
%                 fwhm = [fitresult2G.w1  fitresult2G.w2 ]*2.35;
%                 FF(i,1:2) = [fitresult2G.w1  fitresult2G.w2 ]*2.35;
                             
                h1 = figure(1);
                subplot(3,6,i-1)
                plot(nm,profile,'-o', ...
                    nm,profile_1L,'--r');
                axis([100 900 0 1.05]);
                xlabel(['wR = ' num2str(round(coeffvals_1L(2)))]);
                title([num2str(i-1)])

                param(k,:,z) = [coeffvals_1L]; % [FWHM SNR SBR]coeffvals_LC(3) 1-coeffvals_1L(2)
                k=k+1;
                
%             catch err
%                 disp(strcat(num2str(i),' : No Fit'));
%                 param(k,:,z) = zeros(1,10);
%                 k=k+1;
%             end
        end
    end
    
    




