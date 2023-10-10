function concat_scramble_unconcat_pfamdatabase(inputfile, outputfile)
%
% Last modified on July 31, 2013
% Creating a MSA database with no cognate assumption
% 
% function reference_maker(inputfile1,outputfile1)
%  
% INPUTS: 
%   
%   inputfile1:     MSA of HK
%   intputfile2:    MSA of RR
%   outputfile1:    Scrambled HKRR pairings
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    rng('shuffle');
    dom1 = 1:233;
    dom2 = 234:556;
%     Nhk=64;
%     Nrr=112;
%     [header1, align1] = fastaread(inputfile1);
%     M1=length(align1); 
%     [header2, align2] = fastaread(inputfile2);
% %     M2=length(align2);
%    

    [header, align] = fastaread(inputfile);
    align = char(align);
    align1 = align(:, dom1);
    align2 = align(:, dom2);
    M1 = size(align1,1); 
    M2 = size(align2,1); 

    ordering = randperm(M1);
    
    for i=1:M1
%         %removing warning message from fastawrite for appending file
        warnState = warning; %Save the current warning state
        warning('off','Bioinfo:fastawrite:AppendToFile');

        tempSequence=horzcat(align1(i,:),align2(ordering(i),:));
        
        fastawrite(outputfile,header(i), tempSequence);
        
        
%         fastawrite(outputfile1,header_full(i), align_full{i}(1:Nhk));
%         fastawrite(outputfile2,header_full(i), align_full{i}(Nhk+1:end));
        warning(warnState);
    end
    
%     a=1;
end
