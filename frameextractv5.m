% %%%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& 
% By Shaozhen Song
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Frame = frameextractv5(pixel,coefs,rawname,useref,refname,linenum,framenum)
tic

REF=zeros(pixel,linenum);
Frame = complex(zeros(pixel/2,linenum,framenum,'single'),zeros(pixel/2,linenum,framenum,'single'));
K=1:pixel;KES=polyval(coefs,1:pixel);

    fid=fopen(rawname,'r');
    if useref==0;
        frameCalNum=min(2560,framenum);
        bck = zeros(pixel,linenum);
        for frun=1:frameCalNum
            A = fread(fid,[pixel linenum],'ushort');
            bck = bck + A;
        end
        bck = bck/frameCalNum;
        AMean = mean(bck,2);

        REF=repmat(AMean,1,linenum);
        clear  bck AMean;
    elseif useref==1
        ref = importdata(refname)';
        REF = repmat(ref,1,linenum);
    end
    
    fseek(fid, 0, 'bof');
    fprintf('Total:%d :\\',uint16(framenum));
    percd=0;
    for frun=1:framenum
        percc=uint8(10*frun/framenum);if percc~=percd,fprintf('%d0%%..',percc),percd=percc;end              
            A = fread(fid,[pixel linenum],'ushort');
            SpecC = A - REF;  % Background Subtraction
            clear A;
            SpecC=interp1(K,SpecC,KES,'linenumar',0);  % Interpolation
            Dep=fftshift(fft(SpecC,pixel),1)*2;
            Frame(:,:,frun) = Dep(pixel/2+1:pixel,:);
    end
    fprintf('Done.\n')    
    fclose(fid);
    toc
end