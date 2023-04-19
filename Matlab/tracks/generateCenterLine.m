function center_line = generateCenterLine(cones_blue,cones_yellow)
    b_size = size(cones_blue);
    center_line = zeros(2, b_size(2));
    
    for i = 1:b_size(2)
        cone_yellow_index = getConeYellowIndex(cones_blue(:,i), cones_yellow);
        point = zeros(2,1);
        point(1,1) = (cones_blue(1, i) + cones_yellow(1,cone_yellow_index)) / 2.0;
        point(2,1) = (cones_blue(2, i) + cones_yellow(2,cone_yellow_index)) / 2.0;
    
        center_line(:,i) = point;
    end
end

function cone_yellow_index = getConeYellowIndex(cone_blue, cones_yellow)
    y_size = size(cones_yellow);
    d = inf;
    for i = 1:y_size(2)
        new_d = hypot(cone_blue(1) - cones_yellow(1,i), cone_blue(2) - cones_yellow(2,i));
        if new_d < d
            d = new_d;
            cone_yellow_index = i;
        end
    end
end

