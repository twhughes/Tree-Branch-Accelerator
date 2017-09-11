%read in a command line argument N and plot an order N H-tree
function HTree2(N,size, h, c)
	x = 0; 
    y = 0;
    draw(N, x, y, size, h, N, c);
end

%plot an order n H-tree, centered on (x, y) of the given side length
function draw(n, x, y, size, h, N, c)
	if (n ~= 0)
        drawH(x, y, size, h, n, N, c);
        %compute x- and y-coordinates of the 4 half-size H-trees
        if (n == N)
           x0 = x - size;
           x1 = x + size;
        else
           x0 = x - size*c/2;
           x1 = x + size*c/2;
        end
        y0 = y - size*c/2;
        y1 = y + size*c/2;
        %recursively draw 4 half-size H-trees of order n-1
        draw(n-1, x0, y0, size*c/2, h, N, c); % lower left H-tree
        draw(n-1, x0, y1, size*c/2, h, N, c); % upper left H-tree
        draw(n-1, x1, y0, size*c/2, h, N, c); % lower right H-tree
        draw(n-1, x1, y1, size*c/2, h, N, c); % upper right H-tree
    else
        disp('');
    end
end


function drawH(x, y, size, h, n, N, c)
    % compute the coordinates of the 4 tips of the H
    if (n == N)
       x0 = x - size;
       x1 = x + size;
    else
       x0 = x - size*c/2;
       x1 = x + size*c/2;
    end
	y0 = y - size*c/2;
	y1 = y + size*c/2;
	addBox(x0, x0, y0, y1, h); %/ left vertical segment of the H
	addBox(x1, x1, y0, y1, h); % right vertical segment of the H
	addBox(x0, x1, y,  y,  h); % crossbar segment of the H
end

function addBox(x0, x1, y0, y1, h)
    %takes 4 points of line and widens, returning a box
    if (x0 == x1)
        x0 = x0 - h/2;
        x1 = x1 + h/2;
    elseif (y0 == y1)
        y0 = y0 - h/2;
        y1 = y1 + h/2;
    else
        disp('whoopsies');
    end
    
    X  = min(x0,x1);
    Y  = min(y0,y1);
    XX = max(x0,x1);
    YY = max(y0,y1);
    x0 = X;
    x1 = XX;
    y0 = Y;
    y1 = YY;
    
    b = Box([x0 x1; y0 y1; -1 1]);
    global A;
    A = [A, b];
end

