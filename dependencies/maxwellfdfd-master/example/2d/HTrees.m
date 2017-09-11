function b = getBox(x0,x1,y0,y1,h)
    %takes 4 points of line and widens, returning a box
    if (x0 == x1)
        x0 = x0 - h/2;
        x1 = x1 + h/2;
    elseif (y0 == y1)
        y0 = y0 - h/2;
        y1 = y1 + h/2;
    else
        return;
    end
    b = Box([x0 x1; y0 y1; 0 1]);
end

function draw(boxes, x, y, n, size)
    %takes center of new line and returns 4 points
    if     (mod(n,2) == 0) %even => horizontal
        x0 = x - size/2;
        x1 = x + size/2;
        y0 = y;
        y1 = y;
    elseif (mod(n,2) == 1) %odd  => vertical
        x0 = x;
        x1 = x;
        y0 = y - size/2;
        y1 = y + size/2;
    else
        return;
    end
    boxes = boxes + getBox(x0, x1, y0, y1, h);
end

function boxes = HTree(n, c, size, boxes, x, y)
    draw(boxes, x, y, n, size);
    if (n == 0)
        return;
    end    
    HTree(n-1, c, size*c, boxes, x, y);
end

