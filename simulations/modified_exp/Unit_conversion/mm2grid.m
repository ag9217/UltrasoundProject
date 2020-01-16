function grid_pt =mm2grid(mm_x, mm_y, mm_z)

if mm_x<0
    grid_pt(1)=round((-mm_x+25e-3)*128/50e-3);
    
else 
    grid_pt(1)=round(mm_x*128/50e-3);
end 


if mm_y<0
    grid_pt(2)=round((-mm_y+25e-3)*128/50e-3);
    
else 
    grid_pt(2)=round(mm_y*128/50e-3);
end 


if mm_z<0
    grid_pt(3)=round((-mm_z+25e-3)*128/50e-3);
    
else 
    grid_pt(3)=round(mm_z*128/50e-3);
end 



end 