function DisplayTable(FileName, MatrixNames, resultTab1, resultTab2, AlgName1, AlgName2, isLarge, isTRN, DTIInfo)
% Function for display the main part of LaTeX table codes.
% isLarge and isTRN control the type of the table.
% DTIInfo is needed only for comparing small matrices.
fid = fopen(FileName,'w');
fprintf(fid,'\\centering \\fontsize{8}{10}\\selectfont 	\\begin{tabular}{lcc|ccc|ccc} & &\\multicolumn{1}{l|}{} & \\multicolumn{3}{c|}{%s}  & \\multicolumn{3}{c|}{%s} \\\\ \\hline\n', AlgName1, AlgName2);
if isLarge
    fprintf(fid,'matrix   & $n$  &niter &$\\hat{\\beta}(\\widehat{A}) $&timing &niter &$\\hat{\\beta}(\\widehat{A}) $&timing \\\\ \\hline\n');
else
    fprintf(fid,'matrix   & $n$  & $\\beta(\\widehat{A})$  &niter &correctness &timing &niter &correctness &timing      \\\\ \\hline\n');
end
if isTRN
    if ~strcmp(AlgName2, 'TRN')
        error('When isTRN = true, AlgName1 must be Newton-type method, and AlgName2 must be TRN');
    end
end

matNum = length(MatrixNames);
for count=1:matNum
    % load data and preprocessing
    matName = MatrixNames{count};
    if resultTab1(count,1)~=resultTab2(count,1)
        error('The dimension of two resultTabs are not matched.');
    end
    n = resultTab1(count,1);
    % find exact solution for small problems
    if ~isLarge
        for j=1:length(DTIInfo)
            if strcmp(DTIInfo{j,1},matName)
                DTIexact = DTIInfo{j,2};
                break;
            end
        end
    end
    if ~isLarge
        % In this case, resultTab(:,4) is the success rate; in large problems,
        % it is the calculated DTI.
        if resultTab1(count,4)>=0.75
            correctness1='yes';
        elseif resultTab1(count,4)<=0.25
            correctness1='no';
        else
            correctness1='half';
        end
        
        if resultTab2(count,4)>=0.75
            correctness2='yes';
        elseif resultTab2(count,4)<=0.25
            correctness2='no';
        else
            correctness2='half';
        end
    end
    t1  = resultTab1(count,2);      t2  = resultTab2(count,2);
    if isLarge
        fprintf(fid,'%s   & %d & %d & %.13g  &%.3g & %d & %.13g  &%.3g    ',MatrixNames{count},n,resultTab1(count,3),resultTab1(count,4),t1,resultTab2(count,3),resultTab2(count,4),t2);
    else
        fprintf(fid,'%s   & %d  & %.13g  &%d &%s &%.3g &%d &%s &%.3g   ',MatrixNames{count},n,DTIexact,resultTab1(count,3),correctness1,t1,resultTab2(count,3),correctness2,t2);
    end
    if isTRN
        fprintf(fid,'%.2g',t2/t1);
    end
    fprintf(fid,'\\\\ \\hline\n');
end

fclose(fid);
return;


