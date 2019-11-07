function [k, b] = get_line(x1, y1, x2, y2)
    k = [y2 - y1, -(x2 - x1)];
    b = y2 * (x2 - x1) - x2 * (y2 - y1);   
end

