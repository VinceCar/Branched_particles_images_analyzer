function reference=creating_filter(diameter)

to_build_model=3*diameter+(1-mod(3*diameter,2));
[I,J]=ndgrid(1:to_build_model);
Iv=[to_build_model,ceil(to_build_model/2),1];
Jv=[to_build_model,1,to_build_model];
Image=inpolygon(I,J,Iv,Jv);
xcI=ceil(to_build_model/2);
ycI=ceil(to_build_model/2);
[yr,xr]=find(Image);
if xcI-mean(xr)~=0 || ycI-mean(yr)~=0
    diff_x=mean(xr)-xcI;
    diff_y=mean(yr)-ycI;
    if diff_x>0
    Image=[Image,zeros(size(Image,1),round(diff_x))];
    end
    if diff_y>0
        Image=[Image;zeros(round(diff_y),size(Image,2))];
    end
end
Image=[zeros(size(Image,1),10),Image,zeros(size(Image,1),10)];
Image=[zeros(10,size(Image,2));Image;zeros(10,size(Image,2))];

reference_2=0;
previous=0;
for ang_1=0:5:355
   reference=imrotate(Image,ang_1,'nearest','crop');
   if mod(size(reference,2),2)==0
       reference=[reference,zeros(size(reference,1),1)];
   end
   if mod(size(reference,1),2)==0
       reference=[reference;zeros(1,size(reference,2))];
   end
   [yr,xr]=find(reference);
   xc=round(mean(xr));
   yc=round(mean(yr));

   dif_x=xc-ceil(size(reference,2)/2);
   dif_y=yc-ceil(size(reference,1)/2);
   xr=round(xr+dif_x);
   yr=round(yr+dif_y);
   index=sub2ind(size(reference),yr,xr);
   reference=zeros(size(reference));
   reference(index)=1;
%   reference=bwmorph(reference,'thin',Inf);
    if sum(abs(previous-reference),'all')>0
        reference_2=reference_2+reference;
        previous=reference;
    end
end
[y,x]=find(reference_2);
reference=reference_2(min(y):max(y),min(x):max(x));
end