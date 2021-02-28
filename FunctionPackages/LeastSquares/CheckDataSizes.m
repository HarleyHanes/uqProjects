function data = CheckDataSizes(data)
%CheckDataSizes confirms that data components are column vectors of same
%length
    if size(data.x,2)~=1
        if size(data.x,1)~=1
            error('\n%ix%i matrix inputed for data.x\n',size(data.x))
        else
            warning(['\n%ix%i vector inputed for data.x\n reformating to'...
                '%ix%i vector\n'],size(data.x),size(data.x'))
            data.x=data.x';
        end
    end
    if size(data.y,2)~=1
        if size(data.y,1)~=1
            error('\n%ix%i matrix inputed for data.y\n',size(data.y))
        else
            warning(['\n%ix%i vector inputed for data.y\n reformating to'...
                '%ix%i vector\n'],size(data.y),size(data.y'))
            data.y=data.y';
        end
    end
    if length(data.x)~=length(data.y)
        error('data.x is %ix%i while data.y is %ix%i',size(data.x),size(data.y))
    end
    
end

