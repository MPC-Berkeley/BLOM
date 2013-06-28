function readarray=smap_read_array(id,offset,duration)
%id is uuid
%function returns most current reading
%offset is how much time before now to start data read (in hours)
%duration is duration of data read (in hours)
%S4-16 room temp a71dc6d1-420d-56f9-9b02-f63a5d165fcc
%Cory weather station ec2b82c2-aa68-50ad-8710-12ee8ca63ca7

start_time_str = num2str(1000*ceil(86400 * (now-datenum(0,0,0,offset,0,0)- datenum('01-Jan-1970'))));
end_time_str = num2str(1000*floor(86400 * (now-datenum(0,0,0,offset-duration,0,0)-datenum('01-Jan-1970'))));
data_download = urlread(sprintf('http://new.openbms.org/backend/api/data/uuid/%s?starttime=%s&endtime=%s',id,start_time_str,end_time_str));
parsed_data = parse_json(data_download);
time=zeros(length(parsed_data{1}.Readings),1);
values=zeros(length(parsed_data{1}.Readings),1);
for i=1:length(parsed_data{1}.Readings)
    time(i)=(parsed_data{1}.Readings{i}{1})/1000/86400+datenum('01-Jan-1970');
    values(i)=parsed_data{1}.Readings{i}{2};
end
readarray=timeseries(values,time);
