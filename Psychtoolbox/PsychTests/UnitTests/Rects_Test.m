function success = Rects_Test

success = true;

% create some input
r1 = [0 0 20 20];
r2 = [0 0 40 40];
r3 = [100 200 300 400];
r3a = [100 200 301 401];
point1 = [600 600];
point2 = [150 250];

try
    % test: AdjoinRect
    success = tester(@()AdjoinRect(r1,r3,RectLeft),[80 0 100 20]) && success;
    success = tester(@()AdjoinRect(r1,r3,RectTop),[0 180 20 200]) && success;
    success = tester(@()AdjoinRect(r1,r3,RectRight),[300 0 320 20]) && success;
    success = tester(@()AdjoinRect(r1,r3,RectBottom),[0 400 20 420]) && success;
    success = tester(@()AdjoinRect(r1,r3.',RectBottom),[0 400 20 420]) && success;
    success = tester(@()AdjoinRect(r1.',r3,RectBottom),[0 400 20 420].') && success;
    % vectorized
    success = tester(@()AdjoinRect([r1; r2],r3,RectBottom),[0 400 20 420;0 400 40 440]) && success;
    success = tester(@()AdjoinRect(r1,[r2; r3],RectBottom),[0 40 20 60;0 400 20 420]) && success;
    success = tester(@()AdjoinRect(r1.',[r2; r3],RectBottom),[0 40 20 60;0 400 20 420].') && success;
    success = tester(@()AdjoinRect([r1; r2],[r3; r3a],RectBottom),[0 400 20 420;0 401 40 441]) && success;
    success = tester(@()AdjoinRect([r1; r2],[r3; r3a].',RectBottom),[0 400 20 420;0 401 40 441]) && success;
    
    % test: AlignRect
    success = tester(@()AlignRect(r1,r3,RectLeft),[100 0 120 20]) && success;
    success = tester(@()AlignRect(r1,r3,RectTop),[0 200 20 220]) && success;
    success = tester(@()AlignRect(r1,r3,'right'),[280 0 300 20]) && success;
    success = tester(@()AlignRect(r1,r3,RectBottom),[0 380 20 400]) && success;
    success = tester(@()AlignRect(r1,r3,'center'),[190 290 210 310]) && success;
    success = tester(@()AlignRect(r1,r3a,'centerd'),[190 290 210 310]+.5) && success;
    success = tester(@()AlignRect(r1.',r3a,'centerd'),[190 290 210 310].'+.5) && success;
    success = tester(@()AlignRect(r1,r3a.','centerd'),[190 290 210 310]+.5) && success;
    success = tester(@()AlignRect(r1,r3,RectLeft,RectTop),[100 200 120 220]) && success;
    success = tester(@()AlignRect(r1,r3,RectTop,RectLeft),[100 200 120 220]) && success;
    success = tester(@()AlignRect(r1,r3,RectLeft,RectBottom),[100 380 120 400]) && success;
    success = tester(@()AlignRect(r1.',r3,RectLeft,RectBottom),[100 380 120 400].') && success;
    success = tester(@()AlignRect(r1,r3.',RectLeft,RectBottom),[100 380 120 400]) && success;
    success = tester(@()AlignRect(r1,r3,RectBottom,'left'),[100 380 120 400]) && success;
    success = tester(@()AlignRect(r1,r3,RectRight,RectBottom),[280 380 300 400]) && success;
    success = tester(@()AlignRect(r1,r3,RectBottom,RectRight),[280 380 300 400]) && success;
    success = tester(@()AlignRect(r1,r3,RectRight,'top'),[280 200 300 220]) && success;
    success = tester(@()AlignRect(r1,r3,RectTop,RectRight),[280 200 300 220]) && success;
    success = tester(@()AlignRect(r1,r3,'center',RectLeft),[100 290 120 310]) && success;
    success = tester(@()AlignRect(r1,r3,'center',RectTop),[190 200 210 220]) && success;
    success = tester(@()AlignRect(r1,r3,'center',RectRight),[280 290 300 310]) && success;
    success = tester(@()AlignRect(r1,r3,'center',RectBottom),[190 380 210 400]) && success;
    success = tester(@()AlignRect(r1,r3a,'centerd',RectBottom),[190.5 381 210.5 401]) && success;
    % vectorized
    success = tester(@()AlignRect([r1; r2],r3,'center',RectBottom),[190 380 210 400; 180 360 220 400]) && success;
    success = tester(@()AlignRect([r1; r2].',r3,'center',RectBottom),[190 380 210 400; 180 360 220 400].') && success;
    success = tester(@()AlignRect([r1; r2],r3.','center',RectBottom),[190 380 210 400; 180 360 220 400]) && success;
    success = tester(@()AlignRect(r1,[r2; r3],'center',RectBottom),[10 20 30 40; 190 380 210 400]) && success;
    success = tester(@()AlignRect(r1.',[r2; r3],'center',RectBottom),[10 20 30 40; 190 380 210 400].') && success;
    success = tester(@()AlignRect(r1,[r2; r3].','center',RectBottom),[10 20 30 40; 190 380 210 400]) && success;
    success = tester(@()AlignRect([r1; r2],[r3; r3a],'center',RectBottom),[190 380 210 400; 181 361 221 401]) && success;
    success = tester(@()AlignRect([r1; r2],[r3; r3a],'centerd',RectBottom),[190 380 210 400; 180.5 361 220.5 401]) && success;
    success = tester(@()AlignRect([r1; r2],r3,RectLeft,RectBottom),[100 380 120 400; 100 360 140 400]) && success;
    success = tester(@()AlignRect(r1,[r2; r3],RectLeft,RectBottom),[0 20 20 40; 100 380 120 400]) && success;
    success = tester(@()AlignRect([r1; r2],[r3; r3a],RectLeft,RectBottom),[100 380 120 400; 100 361 140 401]) && success;
    
    % test: CenterRect
    success = tester(@()CenterRect(r1,r3),[190 290 210 310]) && success;
    success = tester(@()CenterRect(r3,r1),[-90 -90 110 110]) && success;
    success = tester(@()CenterRect(r3a,r1),[-91 -91 110 110]) && success;
    success = tester(@()CenterRect(r3a.',r1),[-91 -91 110 110].') && success;
    success = tester(@()CenterRect(r3a,r1.'),[-91 -91 110 110]) && success;
    % vectorized
    success = tester(@()CenterRect([r1; r2],r3),[190 290 210 310; 180 280 220 320]) && success;
    success = tester(@()CenterRect(r1,[r2; r3]),[10 10 30 30; 190 290 210 310]) && success;
    success = tester(@()CenterRect([r1; r2],[r3; r3a]),[190 290 210 310; 181 281 221 321]) && success;
    success = tester(@()CenterRect([r1; r2].',[r3; r3a]),[190 290 210 310; 181 281 221 321].') && success;
    success = tester(@()CenterRect([r1; r2],[r3; r3a].'),[190 290 210 310; 181 281 221 321]) && success;
    % test: CenterRectd
    success = tester(@()CenterRectd(r1,r3a),[190 290 210 310]+.5) && success;
    success = tester(@()CenterRectd(r3a,r1),[-90.5 -90.5 110.5 110.5]) && success;
    success = tester(@()CenterRectd(r3a.',r1),[-90.5 -90.5 110.5 110.5].') && success;
    success = tester(@()CenterRectd(r3a,r1.'),[-90.5 -90.5 110.5 110.5]) && success;
    % vectorized
    success = tester(@()CenterRectd([r1; r2],r3a),[190 290 210 310; 180 280 220 320]+.5) && success;
    success = tester(@()CenterRectd(r1,[r2; r3a]),[10 10 30 30; 190.5 290.5 210.5 310.5]) && success;
    success = tester(@()CenterRectd(r1.',[r2; r3a]),[10 10 30 30; 190.5 290.5 210.5 310.5].') && success;
    success = tester(@()CenterRectd(r1,[r2; r3a].'),[10 10 30 30; 190.5 290.5 210.5 310.5]) && success;
    success = tester(@()CenterRectd([r1; r2],[r3; r3a]),[190 290 210 310; 180.5 280.5 220.5 320.5]) && success;
    
    % test: CenterRectOnPoint
    success = tester(@()CenterRectOnPoint(r1,point1(1),point1(2)),[590 590 610 610]) && success;
    success = tester(@()CenterRectOnPoint(r3,point1(1),point1(2)),[500 500 700 700]) && success;
    success = tester(@()CenterRectOnPoint(r3a,point1(1),point1(2)),[499 499 700 700]) && success;
    success = tester(@()CenterRectOnPoint(r3a.',point1(1),point1(2)),[499 499 700 700].') && success;
    success = tester(@()CenterRectOnPoint(r1+.5,point1(1),point1(2)),[590 590 610 610]-.5) && success;
    success = tester(@()CenterRectOnPoint(r3+.5,point1(1),point1(2)),[500 500 700 700]-.5) && success;
    % vectorized
    success = tester(@()CenterRectOnPoint([r1; r2],point1(1),point1(2)),[590 590 610 610;580 580 620 620]) && success;
    success = tester(@()CenterRectOnPoint(r1,point1(1)+[0 10],point1(2)),[590 590 610 610;600 590 620 610]) && success;
    success = tester(@()CenterRectOnPoint(r1,point1(1),point1(2)+[0 10]),[590 590 610 610;590 600 610 620]) && success;
    success = tester(@()CenterRectOnPoint(r1.',point1(1),point1(2)+[0 10]),[590 590 610 610;590 600 610 620].') && success;
    success = tester(@()CenterRectOnPoint(r1,point1(1),point1(2)+[0 10].'),[590 590 610 610;590 600 610 620]) && success;
    success = tester(@()CenterRectOnPoint([r1; r2],point1(1)+[0 10],point1(2)),[590 590 610 610;590 580 630 620]) && success;
    success = tester(@()CenterRectOnPoint([r1; r2],point1(1),point1(2)+[0 10]),[590 590 610 610;580 590 620 630]) && success;
    success = tester(@()CenterRectOnPoint([r1; r2].',point1(1),point1(2)+[0 10]),[590 590 610 610;580 590 620 630].') && success;
    success = tester(@()CenterRectOnPoint([r1; r2],point1(1)+[0 10],point1(2)+[0 10]),[590 590 610 610;590 590 630 630]) && success;
    % test: CenterRectOnPointd
    success = tester(@()CenterRectOnPointd(r1,point1(1)+.5,point1(2)),[590.5 590 610.5 610]) && success;
    success = tester(@()CenterRectOnPointd(r3,point1(1)+.5,point1(2)),[500.5 500 700.5 700]) && success;
    success = tester(@()CenterRectOnPointd(r3a,point1(1),point1(2)),[499.5 499.5 700.5 700.5]) && success;
    success = tester(@()CenterRectOnPointd(r3a.',point1(1),point1(2)),[499.5 499.5 700.5 700.5].') && success;
    success = tester(@()CenterRectOnPointd(r1+.5,point1(1),point1(2)),[590 590 610 610]) && success;
    success = tester(@()CenterRectOnPointd(r3+.5,point1(1),point1(2)),[500 500 700 700]) && success;
    % vectorized
    success = tester(@()CenterRectOnPointd([r1; r2]+.5,point1(1),point1(2)),[590 590 610 610;580 580 620 620]) && success;
    success = tester(@()CenterRectOnPointd(r1+.5,point1(1)+[0 10],point1(2)),[590 590 610 610;600 590 620 610]) && success;
    success = tester(@()CenterRectOnPointd(r1+.5,point1(1),point1(2)+[0 10]),[590 590 610 610;590 600 610 620]) && success;
    success = tester(@()CenterRectOnPointd(r1.'+.5,point1(1),point1(2)+[0 10]),[590 590 610 610;590 600 610 620].') && success;
    success = tester(@()CenterRectOnPointd(r1+.5,point1(1),point1(2)+[0 10].'),[590 590 610 610;590 600 610 620]) && success;
    success = tester(@()CenterRectOnPointd([r1; r2]+.5,point1(1)+[0 10],point1(2)),[590 590 610 610;590 580 630 620]) && success;
    success = tester(@()CenterRectOnPointd([r1; r2]+.5,point1(1),point1(2)+[0 10]),[590 590 610 610;580 590 620 630]) && success;
    success = tester(@()CenterRectOnPointd([r1; r2].'+.5,point1(1),point1(2)+[0 10]),[590 590 610 610;580 590 620 630].') && success;
    success = tester(@()CenterRectOnPointd([r1; r2]+.5,point1(1)+[0 10]+.5,point1(2)+[0 10]),[590.5 590 610.5 610;590.5 590 630.5 630]) && success;
    
    % test: ClipRect
    success = tester(@()ClipRect(r1,r3),[0 0 0 0]) && success;
    success = tester(@()ClipRect(r1.',r3),[0 0 0 0].') && success;
    success = tester(@()ClipRect(r1,r3.'),[0 0 0 0]) && success;
    success = tester(@()ClipRect(r3,r3a),r3) && success;
    success = tester(@()ClipRect(r3.',r3a),r3.') && success;
    success = tester(@()ClipRect(r3,r3a.'),r3) && success;
    % vectorized
    success = tester(@()ClipRect([r1; r2],r3),[0 0 0 0; 0 0 0 0]) && success;
    success = tester(@()ClipRect(r1,[r3; r3a]),[0 0 0 0; 0 0 0 0]) && success;
    success = tester(@()ClipRect(r1.',[r3; r3a]),[0 0 0 0; 0 0 0 0].') && success;
    success = tester(@()ClipRect(r1,[r3; r3a].'),[0 0 0 0; 0 0 0 0]) && success;
    success = tester(@()ClipRect([r1; r2],[r3; r3a]),[0 0 0 0; 0 0 0 0]) && success;
    success = tester(@()ClipRect([r1; r2].',[r3; r3a]),[0 0 0 0; 0 0 0 0].') && success;
    success = tester(@()ClipRect([r1; r2],[r3; r3a].'),[0 0 0 0; 0 0 0 0]) && success;
    success = tester(@()ClipRect([r1; r3],r3a),[0 0 0 0; r3]) && success;
    success = tester(@()ClipRect(r1,[r2; r3a]),[r1; 0 0 0 0]) && success;
    
    % test: GrowRect
    success = tester(@()GrowRect(r1,0,2),[0 -2 20 22]) && success;
    success = tester(@()GrowRect(r1,2,-2),[-2 2 22 18]) && success;
    success = tester(@()GrowRect(r1.',2,-2),[-2 2 22 18].') && success;
    % vectorized
    success = tester(@()GrowRect(r1,[2 4],-2),[-2 2 22 18;-4 2 24 18]) && success;
    success = tester(@()GrowRect(r1.',[2 4],-2),[-2 2 22 18;-4 2 24 18].') && success;
    success = tester(@()GrowRect([r1; r2],-2,[2 4]),[2 -2 18 22;2 -4 38 44]) && success;
    success = tester(@()GrowRect([r1; r2],-2,4),[2 -4 18 24;2 -4 38 44]) && success;
    success = tester(@()GrowRect([r1; r2].',-2,4),[2 -4 18 24;2 -4 38 44].') && success;
    
    % test: InsetRect
    success = tester(@()InsetRect(r1,0,2),[0 2 20 18]) && success;
    success = tester(@()InsetRect(r1,2,-2),[2 -2 18 22]) && success;
    success = tester(@()InsetRect(r1.',2,-2),[2 -2 18 22].') && success;
    % vectorized
    success = tester(@()InsetRect(r1,[2 4],-2),[2 -2 18 22;4 -2 16 22]) && success;
    success = tester(@()InsetRect(r1.',[2 4],-2),[2 -2 18 22;4 -2 16 22].') && success;
    success = tester(@()InsetRect([r1; r2],-2,[2 4]),[-2 2 22 18;-2 4 42 36]) && success;
    success = tester(@()InsetRect([r1; r2],-2,4),[-2 4 22 16;-2 4 42 36]) && success;
    success = tester(@()InsetRect([r1; r2].',-2,4),[-2 4 22 16;-2 4 42 36].') && success;
    
    % test: IsEmptyRect
    success = tester(@()IsEmptyRect(r1),false) && success;
    success = tester(@()IsEmptyRect(-r1),true) && success;
    success = tester(@()IsEmptyRect(-r1.'),true) && success;
    % vectorized
    success = tester(@()IsEmptyRect([r1; -r1]),[false; true]) && success;
    success = tester(@()IsEmptyRect([r1; -r1].'),[false true]) && success;
    
    % test: IsInRect
    success = tester(@()IsInRect(point1(1),point1(2),r1),false) && success;
    success = tester(@()IsInRect(point2(1),point2(2),r3),true) && success;
    % vectorized
    success = tester(@()IsInRect([point1(1); point2(2)],[point1(2); point2(2)],r3),[false; true]) && success;
    success = tester(@()IsInRect([point1(1); point2(2)].',[point1(2); point2(2)],r3),[false; true]) && success;
    success = tester(@()IsInRect(point2(1),point2(2),[r1;r3]),[false; true]) && success;
    success = tester(@()IsInRect(point2(1),point2(2),[r1;r3].'),[false true]) && success;
    success = tester(@()IsInRect([point1(1); point2(2)],[point1(2); point2(2)],[r3;r3a]),[false; true]) && success;
    success = tester(@()IsInRect([point1(1); point2(2)],[point1(2); point2(2)],[r3;r3a].'),[false true]) && success;
    success = tester(@()IsInRect([point1(1); point2(2)],[point1(2); point2(2)].',[r3;r3a]),[false; true]) && success;
    
    % test: OffsetRect
    success = tester(@()OffsetRect(r1,10,10),[10 10 30 30]) && success;
    success = tester(@()OffsetRect(r1,-10,-20),[-10 -20 10 0]) && success;
    success = tester(@()OffsetRect(r1.',-10,-20),[-10 -20 10 0].') && success;
    % vectorized
    success = tester(@()OffsetRect([r1;r2],10,10),[10 10 30 30;10 10 50 50]) && success;
    success = tester(@()OffsetRect([r1;r2],[10  20],10),[10 10 30 30;20 10 60 50]) && success;
    success = tester(@()OffsetRect([r1;r2],[10; 20],10),[10 10 30 30;20 10 60 50]) && success;
    success = tester(@()OffsetRect([r1;r2].',[10; 20],10),[10 10 30 30;20 10 60 50].') && success;
    success = tester(@()OffsetRect([r1;r2],[10; 20].',10),[10 10 30 30;20 10 60 50]) && success;
    success = tester(@()OffsetRect(r1,[10; 20],10),[10 10 30 30;20 10 40 30]) && success;
    success = tester(@()OffsetRect(r1.',[10; 20],10),[10 10 30 30;20 10 40 30].') && success;
    
    % test: RectCenter
    success = tester(@()RectCenter(r1),10,10) && success;
    success = tester(@()RectCenter(r3a),201,301) && success;
    success = tester(@()RectCenter(r3a.'),201,301) && success;
    % vectorized
    success = tester(@()RectCenter([r1;r2]),[10; 20],[10;20]) && success;
    success = tester(@()RectCenter([r1;r2].'),[10 20],[10 20]) && success;
    
    % test: RectCenter
    success = tester(@()RectCenterd(r1),10,10) && success;
    success = tester(@()RectCenterd(r3a),200.5,300.5) && success;
    success = tester(@()RectCenterd(r3a.'),200.5,300.5) && success;
    % vectorized
    success = tester(@()RectCenterd([r1;r3a]),[10; 200.5],[10;300.5]) && success;
    success = tester(@()RectCenterd([r1;r3a].'),[10 200.5],[10 300.5]) && success;
    
    % test: RectHeight
    success = tester(@()RectHeight(r1),20) && success;
    success = tester(@()RectHeight(r3a),201) && success;
    success = tester(@()RectHeight(r3a.'),201) && success;
    % vectorized
    success = tester(@()RectHeight([r1;r3a]),[20; 201]) && success;
    success = tester(@()RectHeight([r1;r3a].'),[20 201]) && success;
    
    % test: RectOfMatrix
    success = tester(@()RectOfMatrix(zeros(300,100)),[0 0 100 300]) && success;
    success = tester(@()RectOfMatrix(zeros(300,100,3)),[0 0 100 300]) && success;
    
    % test: RectSize
    success = tester(@()RectSize(r1),20,20) && success;
    success = tester(@()RectSize(r3a),201,201) && success;
    success = tester(@()RectSize(r3a.'),201,201) && success;
    % vectorized
    success = tester(@()RectSize([r1;r3a]),[20; 201],[20; 201]) && success;
    success = tester(@()RectSize([r1;r3a].'),[20 201],[20 201]) && success;
    
    % test: RectWidth
    success = tester(@()RectWidth(r1),20) && success;
    success = tester(@()RectWidth(r3a),201) && success;
    success = tester(@()RectWidth(r3a.'),201) && success;
    % vectorized
    success = tester(@()RectWidth([r1;r3a]),[20; 201]) && success;
    success = tester(@()RectWidth([r1;r3a].'),[20 201]) && success;
    
    % test: ScaleRect
    success = tester(@()ScaleRect(r1,1,2),[0 0 20 40]) && success;
    success = tester(@()ScaleRect(r1,2,-2),[0 0 40 -40]) && success;
    success = tester(@()ScaleRect(r1.',2,-2),[0 0 40 -40].') && success;
    % vectorized
    success = tester(@()ScaleRect(r1,[2 4],-2),[0 0 40 -40; 0 0 80 -40]) && success;
    success = tester(@()ScaleRect([r1; r2],-2,[2 4]),[0 0 -40 40; 0 0 -80 160]) && success;
    success = tester(@()ScaleRect([r1; r2].',-2,[2 4]),[0 0 -40 40; 0 0 -80 160].') && success;
    success = tester(@()ScaleRect([r1; r2],-2,4),[0 0 -40 80; 0 0 -80 160]) && success;
    success = tester(@()ScaleRect([r1; r2].',-2,4),[0 0 -40 80; 0 0 -80 160].') && success;
    
    % test: SetRect
    success = tester(@()SetRect(1,2,3,4),[1 2 3 4]) && success;
    success = tester(@()SetRect(1,-2,3,-4),[1 -2 3 -4]) && success;
    % vectorized
    success = tester(@()SetRect([1;2],2,3,4),[1 2 3 4; 2 2 3 4]) && success;
    success = tester(@()SetRect([1 2],2,3,4),[1 2; 2 2; 3 3; 4 4]) && success;
    
    % test: SizeOfRect
    success = tester(@()SizeOfRect(r1),[20 20]) && success;
    success = tester(@()SizeOfRect(r3a),[201 201]) && success;
    success = tester(@()SizeOfRect(r3a.'),[201 201].') && success;
    % vectorized
    success = tester(@()SizeOfRect([r1;r3a]),[20 20; 201 201]) && success;
    success = tester(@()SizeOfRect([r1;r3a].'),[20 201; 20 201]) && success;
    
    % test: UnionRect
    success = tester(@()UnionRect(r1,r3),[0 0 300 400]) && success;
    success = tester(@()UnionRect(r3,r3a),[100 200 301 401]) && success;
    success = tester(@()UnionRect(r3.',r3a),[100 200 301 401].') && success;
    success = tester(@()UnionRect(r3,r3a.'),[100 200 301 401]) && success;
    success = tester(@()UnionRect(r3.',r3a.'),[100 200 301 401].') && success;
    % vectorized
    success = tester(@()UnionRect([r1; r2],r3),[0 0 300 400; 0 0 300 400]) && success;
    success = tester(@()UnionRect(r1,[r3; r3a]),[0 0 300 400; 0 0 301 401]) && success;
    success = tester(@()UnionRect(r1.',[r3; r3a]),[0 0 300 400; 0 0 301 401].') && success;
    success = tester(@()UnionRect(r1,[r3; r3a].'),[0 0 300 400; 0 0 301 401]) && success;
    success = tester(@()UnionRect([r1; r2],[r3; r3a]),[0 0 300 400; 0 0 301 401]) && success;
    success = tester(@()UnionRect([r1; r2].',[r3; r3a]),[0 0 300 400; 0 0 301 401].') && success;
    success = tester(@()UnionRect([r1; r2],[r3; r3a].'),[0 0 300 400; 0 0 301 401]) && success;
    success = tester(@()UnionRect([r1; r3],r3a),[0 0 301 401;100 200 301 401]) && success;
    success = tester(@()UnionRect(r1,[r2; r3a]),[0 0 40 40; 0 0 301 401]) && success;
    
catch me
    success = false;
    fprintf('Unit test %s failed, error ocurred:\n%s\n',mfilename,me.getReport());
end


function [success] = tester(funcall,varargin)

out = cell(size(varargin));
[out{:}] = funcall();

if ~isequal(out,varargin)
    success = false;
    str = func2str(funcall);
    fprintf('%s failed\n',str(4:end));
else
    success = true;
end
