function finaldata = InsertNaN(data, diffthreshold)
% FUNCION  insertNaN: Used to insert NaN values into a vector or matrix
% when the difference value of successive points exceed an input threshold.
%  If vector data is provided, NaN's are inserted where the differences
%  exceed the requested threshold. 
%  If matrix data is input, the difference condition is applied to the
%  first column of data, and NaN's are inserted along the entire row.
%
% usage:  output=insertNaN(data,threshold);
%   INPUTS:     data - input data which will be checked for gaps
%               threshold - threshold value to distinguish where NaN's are
%               placed in the data
%   OUTPUTS:    output - input data, with NaN values inserted where
%               differences exceeded the requested threshold
% 
% Example 1:
% output=insertNaN([11:13 15:17 19:21 25:27],1);
%      returns: 
%      output = [11 12 13 NaN 15 16 17 NaN 19 20 21 NaN 25 26 27]
% Example 2:
%  output=insertNaN([[1:2 5:7 9:10].',[1:7].',[11:17].'],1);
%      returns:
%      output =
%         1     1    11
%         2     2    12
%       NaN   NaN   NaN
%         5     3    13
%         6     4    14
%         7     5    15
%       NaN   NaN   NaN
%         9     6    16
%        10     7    17
% 
% Chris Miller
%cwmiller@nps.edu
% 9/14/11
if isvector(data), 
    diffdata=diff(data);
    index=find(diff(data)>diffthreshold);
    if isempty(index),
        finaldata=data;
        return;
    end;
    finaldata=NaN*ones(1,length(data)+length(index));  % preallocate output 
    finaldata(1:index(1))=data(1:index(1));
    if length(index)>1,
        for i=2:length(index),
            finaldata(index(i-1)+i:index(i)+i-1)=data(index(i-1)+1:index(i));
        end;
    else
        i=1;
    end;
    finaldata(index(i)+i+1:length(finaldata))=data(index(i)+1:length(data));
else,
    diffdata=diff(data(:,1));
    index=find(diffdata>diffthreshold);
    if isempty(index),
        finaldata=data;
        return;
    end;
    [n,m]=size(data);
    finaldata=NaN*ones(n+length(index),m);  % preallocate output 
    finaldata(1:index(1),:)=data(1:index(1),:);
    if length(index)>1,
        for i=2:length(index),
            finaldata(index(i-1)+i:index(i)+i-1,:)=data(index(i-1)+1:index(i),:);
        end;
    else
        i=1;
    end;
    finaldata(index(i)+i+1:length(finaldata),:)=data(index(i)+1:length(data),:);
end;