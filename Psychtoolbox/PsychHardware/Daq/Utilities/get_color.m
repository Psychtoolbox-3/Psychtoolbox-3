function colorvector = get_color(colorstring)
% Syntax: ColorVector = get_color(ColorNameString)
%
% convert color name to color triple. Color name/triple combinations
% were taken from the Unix X11/R4 rgb.txt file and converted to
% a format matlab would properly interpret.  This function was
% created by Mickey P. Rowe some time around 1993 or 1994.  In June of 2000
% I found out that the switch statement could be used instead of successive 
% calls to strcmp.  I suspect that the compiled version is identical, but this
% looks more elegant.

while colorstring(length(colorstring)) == setstr(' '),
  colorstring(length(colorstring)) = [];
end;

switch(colorstring)
	case 'snow',
	 colorvector = [ 1 0.98039216 0.98039216 ];
	case 'ghost white',
	 colorvector = [ 0.97254902 0.97254902 1 ];
	case 'GhostWhite',
	 colorvector = [ 0.97254902 0.97254902 1 ];
	case 'white smoke',
	 colorvector = [ 0.96078431 0.96078431 0.96078431 ];
	case 'WhiteSmoke',
	 colorvector = [ 0.96078431 0.96078431 0.96078431 ];
	case 'gainsboro',
	 colorvector = [ 0.86274510 0.86274510 0.86274510 ];
	case 'floral white',
	 colorvector = [ 1 0.98039216 0.94117647 ];
	case 'FloralWhite',
	 colorvector = [ 1 0.98039216 0.94117647 ];
	case 'old lace',
	 colorvector = [ 0.99215686 0.96078431 0.90196078 ];
	case 'OldLace',
	 colorvector = [ 0.99215686 0.96078431 0.90196078 ];
	case 'linen',
	 colorvector = [ 0.98039216 0.94117647 0.90196078 ];
	case 'antique white',
	 colorvector = [ 0.98039216 0.92156863 0.84313725 ];
	case 'AntiqueWhite',
	 colorvector = [ 0.98039216 0.92156863 0.84313725 ];
	case 'papaya whip',
	 colorvector = [ 1 0.93725490 0.83529412 ];
	case 'PapayaWhip',
	 colorvector = [ 1 0.93725490 0.83529412 ];
	case 'blanched almond',
	 colorvector = [ 1 0.92156863 0.80392157 ];
	case 'BlanchedAlmond',
	 colorvector = [ 1 0.92156863 0.80392157 ];
	case 'bisque',
	 colorvector = [ 1 0.89411765 0.76862745 ];
	case 'peach puff',
	 colorvector = [ 1 0.85490196 0.72549020 ];
	case 'PeachPuff',
	 colorvector = [ 1 0.85490196 0.72549020 ];
	case 'navajo white',
	 colorvector = [ 1 0.87058824 0.67843137 ];
	case 'NavajoWhite',
	 colorvector = [ 1 0.87058824 0.67843137 ];
	case 'moccasin',
	 colorvector = [ 1 0.89411765 0.70980392 ];
	case 'cornsilk',
	 colorvector = [ 1 0.97254902 0.86274510 ];
	case 'ivory',
	 colorvector = [ 1 1 0.94117647 ];
	case 'lemon chiffon',
	 colorvector = [ 1 0.98039216 0.80392157 ];
	case 'LemonChiffon',
	 colorvector = [ 1 0.98039216 0.80392157 ];
	case 'seashell',
	 colorvector = [ 1 0.96078431 0.93333333 ];
	case 'honeydew',
	 colorvector = [ 0.94117647 1 0.94117647 ];
	case 'mint cream',
	 colorvector = [ 0.96078431 1 0.98039216 ];
	case 'MintCream',
	 colorvector = [ 0.96078431 1 0.98039216 ];
	case 'azure',
	 colorvector = [ 0.94117647 1 1 ];
	case 'alice blue',
	 colorvector = [ 0.94117647 0.97254902 1 ];
	case 'AliceBlue',
	 colorvector = [ 0.94117647 0.97254902 1 ];
	case 'lavender',
	 colorvector = [ 0.90196078 0.90196078 0.98039216 ];
	case 'lavender blush',
	 colorvector = [ 1 0.94117647 0.96078431 ];
	case 'LavenderBlush',
	 colorvector = [ 1 0.94117647 0.96078431 ];
	case 'misty rose',
	 colorvector = [ 1 0.89411765 0.88235294 ];
	case 'MistyRose',
	 colorvector = [ 1 0.89411765 0.88235294 ];
	case 'white',
	 colorvector = [ 1 1 1 ];
	case 'black',
	 colorvector = [ 0 0 0 ];
	case 'dark slate gray',
	 colorvector = [ 0.18431373 0.30980392 0.30980392 ];
	case 'DarkSlateGray',
	 colorvector = [ 0.18431373 0.30980392 0.30980392 ];
	case 'dark slate grey',
	 colorvector = [ 0.18431373 0.30980392 0.30980392 ];
	case 'DarkSlateGrey',
	 colorvector = [ 0.18431373 0.30980392 0.30980392 ];
	case 'dim gray',
	 colorvector = [ 0.41176471 0.41176471 0.41176471 ];
	case 'DimGray',
	 colorvector = [ 0.41176471 0.41176471 0.41176471 ];
	case 'dim grey',
	 colorvector = [ 0.41176471 0.41176471 0.41176471 ];
	case 'DimGrey',
	 colorvector = [ 0.41176471 0.41176471 0.41176471 ];
	case 'slate gray',
	 colorvector = [ 0.43921569 0.50196078 0.56470588 ];
	case 'SlateGray',
	 colorvector = [ 0.43921569 0.50196078 0.56470588 ];
	case 'slate grey',
	 colorvector = [ 0.43921569 0.50196078 0.56470588 ];
	case 'SlateGrey',
	 colorvector = [ 0.43921569 0.50196078 0.56470588 ];
	case 'light slate gray',
	 colorvector = [ 0.46666667 0.53333333 0.6 ];
	case 'LightSlateGray',
	 colorvector = [ 0.46666667 0.53333333 0.6 ];
	case 'light slate grey',
	 colorvector = [ 0.46666667 0.53333333 0.6 ];
	case 'LightSlateGrey',
	 colorvector = [ 0.46666667 0.53333333 0.6 ];
	case 'gray',
	 colorvector = [ 0.74509804 0.74509804 0.74509804 ];
	case 'grey',
	 colorvector = [ 0.74509804 0.74509804 0.74509804 ];
	case 'light grey',
	 colorvector = [ 0.82745098 0.82745098 0.82745098 ];
	case 'LightGrey',
	 colorvector = [ 0.82745098 0.82745098 0.82745098 ];
	case 'light gray',
	 colorvector = [ 0.82745098 0.82745098 0.82745098 ];
	case 'LightGray',
	 colorvector = [ 0.82745098 0.82745098 0.82745098 ];
	case 'midnight blue',
	 colorvector = [ 0.098039216 0.098039216 0.43921569 ];
	case 'MidnightBlue',
	 colorvector = [ 0.098039216 0.098039216 0.43921569 ];
	case 'navy',
	 colorvector = [ 0 0 0.50196078 ];
	case 'navy blue',
	 colorvector = [ 0 0 0.50196078 ];
	case 'NavyBlue',
	 colorvector = [ 0 0 0.50196078 ];
	case 'cornflower blue',
	 colorvector = [ 0.39215686 0.58431373 0.92941176 ];
	case 'CornflowerBlue',
	 colorvector = [ 0.39215686 0.58431373 0.92941176 ];
	case 'dark slate blue',
	 colorvector = [ 0.28235294 0.23921569 0.54509804 ];
	case 'DarkSlateBlue',
	 colorvector = [ 0.28235294 0.23921569 0.54509804 ];
	case 'slate blue',
	 colorvector = [ 0.41568627 0.35294118 0.80392157 ];
	case 'SlateBlue',
	 colorvector = [ 0.41568627 0.35294118 0.80392157 ];
	case 'medium slate blue',
	 colorvector = [ 0.48235294 0.40784314 0.93333333 ];
	case 'MediumSlateBlue',
	 colorvector = [ 0.48235294 0.40784314 0.93333333 ];
	case 'light slate blue',
	 colorvector = [ 0.51764706 0.43921569 1 ];
	case 'LightSlateBlue',
	 colorvector = [ 0.51764706 0.43921569 1 ];
	case 'medium blue',
	 colorvector = [ 0 0 0.80392157 ];
	case 'MediumBlue',
	 colorvector = [ 0 0 0.80392157 ];
	case 'royal blue',
	 colorvector = [ 0.25490196 0.41176471 0.88235294 ];
	case 'RoyalBlue',
	 colorvector = [ 0.25490196 0.41176471 0.88235294 ];
	case 'blue',
	 colorvector = [ 0 0 1 ];
	case 'dodger blue',
	 colorvector = [ 0.11764706 0.56470588 1 ];
	case 'DodgerBlue',
	 colorvector = [ 0.11764706 0.56470588 1 ];
	case 'deep sky blue',
	 colorvector = [ 0 0.74901961 1 ];
	case 'DeepSkyBlue',
	 colorvector = [ 0 0.74901961 1 ];
	case 'sky blue',
	 colorvector = [ 0.52941176 0.80784314 0.92156863 ];
	case 'SkyBlue',
	 colorvector = [ 0.52941176 0.80784314 0.92156863 ];
	case 'light sky blue',
	 colorvector = [ 0.52941176 0.80784314 0.98039216 ];
	case 'LightSkyBlue',
	 colorvector = [ 0.52941176 0.80784314 0.98039216 ];
	case 'steel blue',
	 colorvector = [ 0.27450980 0.50980392 0.70588235 ];
	case 'SteelBlue',
	 colorvector = [ 0.27450980 0.50980392 0.70588235 ];
	case 'light steel blue',
	 colorvector = [ 0.69019608 0.76862745 0.87058824 ];
	case 'LightSteelBlue',
	 colorvector = [ 0.69019608 0.76862745 0.87058824 ];
	case 'light blue',
	 colorvector = [ 0.67843137 0.84705882 0.90196078 ];
	case 'LightBlue',
	 colorvector = [ 0.67843137 0.84705882 0.90196078 ];
	case 'powder blue',
	 colorvector = [ 0.69019608 0.87843137 0.90196078 ];
	case 'PowderBlue',
	 colorvector = [ 0.69019608 0.87843137 0.90196078 ];
	case 'pale turquoise',
	 colorvector = [ 0.68627451 0.93333333 0.93333333 ];
	case 'PaleTurquoise',
	 colorvector = [ 0.68627451 0.93333333 0.93333333 ];
	case 'dark turquoise',
	 colorvector = [ 0 0.80784314 0.81960784 ];
	case 'DarkTurquoise',
	 colorvector = [ 0 0.80784314 0.81960784 ];
	case 'medium turquoise',
	 colorvector = [ 0.28235294 0.81960784 0.8 ];
	case 'MediumTurquoise',
	 colorvector = [ 0.28235294 0.81960784 0.8 ];
	case 'turquoise',
	 colorvector = [ 0.25098039 0.87843137 0.81568627 ];
	case 'cyan',
	 colorvector = [ 0 1 1 ];
	case 'light cyan',
	 colorvector = [ 0.87843137 1 1 ];
	case 'LightCyan',
	 colorvector = [ 0.87843137 1 1 ];
	case 'cadet blue',
	 colorvector = [ 0.37254902 0.61960784 0.62745098 ];
	case 'CadetBlue',
	 colorvector = [ 0.37254902 0.61960784 0.62745098 ];
	case 'medium aquamarine',
	 colorvector = [ 0.4 0.80392157 0.66666667 ];
	case 'MediumAquamarine',
	 colorvector = [ 0.4 0.80392157 0.66666667 ];
	case 'aquamarine',
	 colorvector = [ 0.49803922 1 0.83137255 ];
	case 'dark green',
	 colorvector = [ 0 0.39215686 0 ];
	case 'DarkGreen',
	 colorvector = [ 0 0.39215686 0 ];
	case 'dark olive green',
	 colorvector = [ 0.33333333 0.41960784 0.18431373 ];
	case 'DarkOliveGreen',
	 colorvector = [ 0.33333333 0.41960784 0.18431373 ];
	case 'dark sea green',
	 colorvector = [ 0.56078431 0.73725490 0.56078431 ];
	case 'DarkSeaGreen',
	 colorvector = [ 0.56078431 0.73725490 0.56078431 ];
	case 'sea green',
	 colorvector = [ 0.18039216 0.54509804 0.34117647 ];
	case 'SeaGreen',
	 colorvector = [ 0.18039216 0.54509804 0.34117647 ];
	case 'medium sea green',
	 colorvector = [ 0.23529412 0.70196078 0.44313725 ];
	case 'MediumSeaGreen',
	 colorvector = [ 0.23529412 0.70196078 0.44313725 ];
	case 'light sea green',
	 colorvector = [ 0.12549020 0.69803922 0.66666667 ];
	case 'LightSeaGreen',
	 colorvector = [ 0.12549020 0.69803922 0.66666667 ];
	case 'pale green',
	 colorvector = [ 0.59607843 0.98431373 0.59607843 ];
	case 'PaleGreen',
	 colorvector = [ 0.59607843 0.98431373 0.59607843 ];
	case 'spring green',
	 colorvector = [ 0 1 0.49803922 ];
	case 'SpringGreen',
	 colorvector = [ 0 1 0.49803922 ];
	case 'lawn green',
	 colorvector = [ 0.48627451 0.98823529 0 ];
	case 'LawnGreen',
	 colorvector = [ 0.48627451 0.98823529 0 ];
	case 'green',
	 colorvector = [ 0 1 0 ];
	case 'chartreuse',
	 colorvector = [ 0.49803922 1 0 ];
	case 'medium spring green',
	 colorvector = [ 0 0.98039216 0.60392157 ];
	case 'MediumSpringGreen',
	 colorvector = [ 0 0.98039216 0.60392157 ];
	case 'green yellow',
	 colorvector = [ 0.67843137 1 0.18431373 ];
	case 'GreenYellow',
	 colorvector = [ 0.67843137 1 0.18431373 ];
	case 'lime green',
	 colorvector = [ 0.19607843 0.80392157 0.19607843 ];
	case 'LimeGreen',
	 colorvector = [ 0.19607843 0.80392157 0.19607843 ];
	case 'yellow green',
	 colorvector = [ 0.60392157 0.80392157 0.19607843 ];
	case 'YellowGreen',
	 colorvector = [ 0.60392157 0.80392157 0.19607843 ];
	case 'forest green',
	 colorvector = [ 0.13333333 0.54509804 0.13333333 ];
	case 'ForestGreen',
	 colorvector = [ 0.13333333 0.54509804 0.13333333 ];
	case 'olive drab',
	 colorvector = [ 0.41960784 0.55686275 0.13725490 ];
	case 'OliveDrab',
	 colorvector = [ 0.41960784 0.55686275 0.13725490 ];
	case 'dark khaki',
	 colorvector = [ 0.74117647 0.71764706 0.41960784 ];
	case 'DarkKhaki',
	 colorvector = [ 0.74117647 0.71764706 0.41960784 ];
	case 'khaki',
	 colorvector = [ 0.94117647 0.90196078 0.54901961 ];
	case 'pale goldenrod',
	 colorvector = [ 0.93333333 0.90980392 0.66666667 ];
	case 'PaleGoldenrod',
	 colorvector = [ 0.93333333 0.90980392 0.66666667 ];
	case 'light goldenrod yellow',
	 colorvector = [ 0.98039216 0.98039216 0.82352941 ];
	case 'LightGoldenrodYellow',
	 colorvector = [ 0.98039216 0.98039216 0.82352941 ];
	case 'light yellow',
	 colorvector = [ 1 1 0.87843137 ];
	case 'LightYellow',
	 colorvector = [ 1 1 0.87843137 ];
	case 'yellow',
	 colorvector = [ 1 1 0 ];
	case 'gold',
	 colorvector = [ 1 0.84313725 0 ];
	case 'light goldenrod',
	 colorvector = [ 0.93333333 0.86666667 0.50980392 ];
	case 'LightGoldenrod',
	 colorvector = [ 0.93333333 0.86666667 0.50980392 ];
	case 'goldenrod',
	 colorvector = [ 0.85490196 0.64705882 0.12549020 ];
	case 'dark goldenrod',
	 colorvector = [ 0.72156863 0.52549020 0.043137255 ];
	case 'DarkGoldenrod',
	 colorvector = [ 0.72156863 0.52549020 0.043137255 ];
	case 'rosy brown',
	 colorvector = [ 0.73725490 0.56078431 0.56078431 ];
	case 'RosyBrown',
	 colorvector = [ 0.73725490 0.56078431 0.56078431 ];
	case 'indian red',
	 colorvector = [ 0.80392157 0.36078431 0.36078431 ];
	case 'IndianRed',
	 colorvector = [ 0.80392157 0.36078431 0.36078431 ];
	case 'saddle brown',
	 colorvector = [ 0.54509804 0.27058824 0.074509804 ];
	case 'SaddleBrown',
	 colorvector = [ 0.54509804 0.27058824 0.074509804 ];
	case 'sienna',
	 colorvector = [ 0.62745098 0.32156863 0.17647059 ];
	case 'peru',
	 colorvector = [ 0.80392157 0.52156863 0.24705882 ];
	case 'burlywood',
	 colorvector = [ 0.87058824 0.72156863 0.52941176 ];
	case 'beige',
	 colorvector = [ 0.96078431 0.96078431 0.86274510 ];
	case 'wheat',
	 colorvector = [ 0.96078431 0.87058824 0.70196078 ];
	case 'sandy brown',
	 colorvector = [ 0.95686275 0.64313725 0.37647059 ];
	case 'SandyBrown',
	 colorvector = [ 0.95686275 0.64313725 0.37647059 ];
	case 'tan',
	 colorvector = [ 0.82352941 0.70588235 0.54901961 ];
	case 'chocolate',
	 colorvector = [ 0.82352941 0.41176471 0.11764706 ];
	case 'firebrick',
	 colorvector = [ 0.69803922 0.13333333 0.13333333 ];
	case 'brown',
	 colorvector = [ 0.64705882 0.16470588 0.16470588 ];
	case 'dark salmon',
	 colorvector = [ 0.91372549 0.58823529 0.47843137 ];
	case 'DarkSalmon',
	 colorvector = [ 0.91372549 0.58823529 0.47843137 ];
	case 'salmon',
	 colorvector = [ 0.98039216 0.50196078 0.44705882 ];
	case 'light salmon',
	 colorvector = [ 1 0.62745098 0.47843137 ];
	case 'LightSalmon',
	 colorvector = [ 1 0.62745098 0.47843137 ];
	case 'orange',
	 colorvector = [ 1 0.64705882 0 ];
	case 'dark orange',
	 colorvector = [ 1 0.54901961 0 ];
	case 'DarkOrange',
	 colorvector = [ 1 0.54901961 0 ];
	case 'coral',
	 colorvector = [ 1 0.49803922 0.31372549 ];
	case 'light coral',
	 colorvector = [ 0.94117647 0.50196078 0.50196078 ];
	case 'LightCoral',
	 colorvector = [ 0.94117647 0.50196078 0.50196078 ];
	case 'tomato',
	 colorvector = [ 1 0.38823529 0.27843137 ];
	case 'orange red',
	 colorvector = [ 1 0.27058824 0 ];
	case 'OrangeRed',
	 colorvector = [ 1 0.27058824 0 ];
	case 'red',
	 colorvector = [ 1 0 0 ];
	case 'hot pink',
	 colorvector = [ 1 0.41176471 0.70588235 ];
	case 'HotPink',
	 colorvector = [ 1 0.41176471 0.70588235 ];
	case 'deep pink',
	 colorvector = [ 1 0.078431373 0.57647059 ];
	case 'DeepPink',
	 colorvector = [ 1 0.078431373 0.57647059 ];
	case 'pink',
	 colorvector = [ 1 0.75294118 0.79607843 ];
	case 'light pink',
	 colorvector = [ 1 0.71372549 0.75686275 ];
	case 'LightPink',
	 colorvector = [ 1 0.71372549 0.75686275 ];
	case 'pale violet red',
	 colorvector = [ 0.85882353 0.43921569 0.57647059 ];
	case 'PaleVioletRed',
	 colorvector = [ 0.85882353 0.43921569 0.57647059 ];
	case 'maroon',
	 colorvector = [ 0.69019608 0.18823529 0.37647059 ];
	case 'medium violet red',
	 colorvector = [ 0.78039216 0.082352941 0.52156863 ];
	case 'MediumVioletRed',
	 colorvector = [ 0.78039216 0.082352941 0.52156863 ];
	case 'violet red',
	 colorvector = [ 0.81568627 0.12549020 0.56470588 ];
	case 'VioletRed',
	 colorvector = [ 0.81568627 0.12549020 0.56470588 ];
	case 'magenta',
	 colorvector = [ 1 0 1 ];
	case 'violet',
	 colorvector = [ 0.93333333 0.50980392 0.93333333 ];
	case 'plum',
	 colorvector = [ 0.86666667 0.62745098 0.86666667 ];
	case 'orchid',
	 colorvector = [ 0.85490196 0.43921569 0.83921569 ];
	case 'medium orchid',
	 colorvector = [ 0.72941176 0.33333333 0.82745098 ];
	case 'MediumOrchid',
	 colorvector = [ 0.72941176 0.33333333 0.82745098 ];
	case 'dark orchid',
	 colorvector = [ 0.6 0.19607843 0.8 ];
	case 'DarkOrchid',
	 colorvector = [ 0.6 0.19607843 0.8 ];
	case 'dark violet',
	 colorvector = [ 0.58039216 0 0.82745098 ];
	case 'DarkViolet',
	 colorvector = [ 0.58039216 0 0.82745098 ];
	case 'blue violet',
	 colorvector = [ 0.54117647 0.16862745 0.88627451 ];
	case 'BlueViolet',
	 colorvector = [ 0.54117647 0.16862745 0.88627451 ];
	case 'purple',
	 colorvector = [ 0.62745098 0.12549020 0.94117647 ];
	case 'medium purple',
	 colorvector = [ 0.57647059 0.43921569 0.85882353 ];
	case 'MediumPurple',
	 colorvector = [ 0.57647059 0.43921569 0.85882353 ];
	case 'thistle',
	 colorvector = [ 0.84705882 0.74901961 0.84705882 ];
	case 'snow1',
	 colorvector = [ 1 0.98039216 0.98039216 ];
	case 'snow2',
	 colorvector = [ 0.93333333 0.91372549 0.91372549 ];
	case 'snow3',
	 colorvector = [ 0.80392157 0.78823529 0.78823529 ];
	case 'snow4',
	 colorvector = [ 0.54509804 0.53725490 0.53725490 ];
	case 'seashell1',
	 colorvector = [ 1 0.96078431 0.93333333 ];
	case 'seashell2',
	 colorvector = [ 0.93333333 0.89803922 0.87058824 ];
	case 'seashell3',
	 colorvector = [ 0.80392157 0.77254902 0.74901961 ];
	case 'seashell4',
	 colorvector = [ 0.54509804 0.52549020 0.50980392 ];
	case 'AntiqueWhite1',
	 colorvector = [ 1 0.93725490 0.85882353 ];
	case 'AntiqueWhite2',
	 colorvector = [ 0.93333333 0.87450980 0.8 ];
	case 'AntiqueWhite3',
	 colorvector = [ 0.80392157 0.75294118 0.69019608 ];
	case 'AntiqueWhite4',
	 colorvector = [ 0.54509804 0.51372549 0.47058824 ];
	case 'bisque1',
	 colorvector = [ 1 0.89411765 0.76862745 ];
	case 'bisque2',
	 colorvector = [ 0.93333333 0.83529412 0.71764706 ];
	case 'bisque3',
	 colorvector = [ 0.80392157 0.71764706 0.61960784 ];
	case 'bisque4',
	 colorvector = [ 0.54509804 0.49019608 0.41960784 ];
	case 'PeachPuff1',
	 colorvector = [ 1 0.85490196 0.72549020 ];
	case 'PeachPuff2',
	 colorvector = [ 0.93333333 0.79607843 0.67843137 ];
	case 'PeachPuff3',
	 colorvector = [ 0.80392157 0.68627451 0.58431373 ];
	case 'PeachPuff4',
	 colorvector = [ 0.54509804 0.46666667 0.39607843 ];
	case 'NavajoWhite1',
	 colorvector = [ 1 0.87058824 0.67843137 ];
	case 'NavajoWhite2',
	 colorvector = [ 0.93333333 0.81176471 0.63137255 ];
	case 'NavajoWhite3',
	 colorvector = [ 0.80392157 0.70196078 0.54509804 ];
	case 'NavajoWhite4',
	 colorvector = [ 0.54509804 0.47450980 0.36862745 ];
	case 'LemonChiffon1',
	 colorvector = [ 1 0.98039216 0.80392157 ];
	case 'LemonChiffon2',
	 colorvector = [ 0.93333333 0.91372549 0.74901961 ];
	case 'LemonChiffon3',
	 colorvector = [ 0.80392157 0.78823529 0.64705882 ];
	case 'LemonChiffon4',
	 colorvector = [ 0.54509804 0.53725490 0.43921569 ];
	case 'cornsilk1',
	 colorvector = [ 1 0.97254902 0.86274510 ];
	case 'cornsilk2',
	 colorvector = [ 0.93333333 0.90980392 0.80392157 ];
	case 'cornsilk3',
	 colorvector = [ 0.80392157 0.78431373 0.69411765 ];
	case 'cornsilk4',
	 colorvector = [ 0.54509804 0.53333333 0.47058824 ];
	case 'ivory1',
	 colorvector = [ 1 1 0.94117647 ];
	case 'ivory2',
	 colorvector = [ 0.93333333 0.93333333 0.87843137 ];
	case 'ivory3',
	 colorvector = [ 0.80392157 0.80392157 0.75686275 ];
	case 'ivory4',
	 colorvector = [ 0.54509804 0.54509804 0.51372549 ];
	case 'honeydew1',
	 colorvector = [ 0.94117647 1 0.94117647 ];
	case 'honeydew2',
	 colorvector = [ 0.87843137 0.93333333 0.87843137 ];
	case 'honeydew3',
	 colorvector = [ 0.75686275 0.80392157 0.75686275 ];
	case 'honeydew4',
	 colorvector = [ 0.51372549 0.54509804 0.51372549 ];
	case 'LavenderBlush1',
	 colorvector = [ 1 0.94117647 0.96078431 ];
	case 'LavenderBlush2',
	 colorvector = [ 0.93333333 0.87843137 0.89803922 ];
	case 'LavenderBlush3',
	 colorvector = [ 0.80392157 0.75686275 0.77254902 ];
	case 'LavenderBlush4',
	 colorvector = [ 0.54509804 0.51372549 0.52549020 ];
	case 'MistyRose1',
	 colorvector = [ 1 0.89411765 0.88235294 ];
	case 'MistyRose2',
	 colorvector = [ 0.93333333 0.83529412 0.82352941 ];
	case 'MistyRose3',
	 colorvector = [ 0.80392157 0.71764706 0.70980392 ];
	case 'MistyRose4',
	 colorvector = [ 0.54509804 0.49019608 0.48235294 ];
	case 'azure1',
	 colorvector = [ 0.94117647 1 1 ];
	case 'azure2',
	 colorvector = [ 0.87843137 0.93333333 0.93333333 ];
	case 'azure3',
	 colorvector = [ 0.75686275 0.80392157 0.80392157 ];
	case 'azure4',
	 colorvector = [ 0.51372549 0.54509804 0.54509804 ];
	case 'SlateBlue1',
	 colorvector = [ 0.51372549 0.43529412 1 ];
	case 'SlateBlue2',
	 colorvector = [ 0.47843137 0.40392157 0.93333333 ];
	case 'SlateBlue3',
	 colorvector = [ 0.41176471 0.34901961 0.80392157 ];
	case 'SlateBlue4',
	 colorvector = [ 0.27843137 0.23529412 0.54509804 ];
	case 'RoyalBlue1',
	 colorvector = [ 0.28235294 0.46274510 1 ];
	case 'RoyalBlue2',
	 colorvector = [ 0.26274510 0.43137255 0.93333333 ];
	case 'RoyalBlue3',
	 colorvector = [ 0.22745098 0.37254902 0.80392157 ];
	case 'RoyalBlue4',
	 colorvector = [ 0.15294118 0.25098039 0.54509804 ];
	case 'blue1',
	 colorvector = [ 0 0 1 ];
	case 'blue2',
	 colorvector = [ 0 0 0.93333333 ];
	case 'blue3',
	 colorvector = [ 0 0 0.80392157 ];
	case 'blue4',
	 colorvector = [ 0 0 0.54509804 ];
	case 'DodgerBlue1',
	 colorvector = [ 0.11764706 0.56470588 1 ];
	case 'DodgerBlue2',
	 colorvector = [ 0.10980392 0.52549020 0.93333333 ];
	case 'DodgerBlue3',
	 colorvector = [ 0.094117647 0.45490196 0.80392157 ];
	case 'DodgerBlue4',
	 colorvector = [ 0.062745098 0.30588235 0.54509804 ];
	case 'SteelBlue1',
	 colorvector = [ 0.38823529 0.72156863 1 ];
	case 'SteelBlue2',
	 colorvector = [ 0.36078431 0.67450980 0.93333333 ];
	case 'SteelBlue3',
	 colorvector = [ 0.30980392 0.58039216 0.80392157 ];
	case 'SteelBlue4',
	 colorvector = [ 0.21176471 0.39215686 0.54509804 ];
	case 'DeepSkyBlue1',
	 colorvector = [ 0 0.74901961 1 ];
	case 'DeepSkyBlue2',
	 colorvector = [ 0 0.69803922 0.93333333 ];
	case 'DeepSkyBlue3',
	 colorvector = [ 0 0.60392157 0.80392157 ];
	case 'DeepSkyBlue4',
	 colorvector = [ 0 0.40784314 0.54509804 ];
	case 'SkyBlue1',
	 colorvector = [ 0.52941176 0.80784314 1 ];
	case 'SkyBlue2',
	 colorvector = [ 0.49411765 0.75294118 0.93333333 ];
	case 'SkyBlue3',
	 colorvector = [ 0.42352941 0.65098039 0.80392157 ];
	case 'SkyBlue4',
	 colorvector = [ 0.29019608 0.43921569 0.54509804 ];
	case 'LightSkyBlue1',
	 colorvector = [ 0.69019608 0.88627451 1 ];
	case 'LightSkyBlue2',
	 colorvector = [ 0.64313725 0.82745098 0.93333333 ];
	case 'LightSkyBlue3',
	 colorvector = [ 0.55294118 0.71372549 0.80392157 ];
	case 'LightSkyBlue4',
	 colorvector = [ 0.37647059 0.48235294 0.54509804 ];
	case 'SlateGray1',
	 colorvector = [ 0.77647059 0.88627451 1 ];
	case 'SlateGray2',
	 colorvector = [ 0.72549020 0.82745098 0.93333333 ];
	case 'SlateGray3',
	 colorvector = [ 0.62352941 0.71372549 0.80392157 ];
	case 'SlateGray4',
	 colorvector = [ 0.42352941 0.48235294 0.54509804 ];
	case 'LightSteelBlue1',
	 colorvector = [ 0.79215686 0.88235294 1 ];
	case 'LightSteelBlue2',
	 colorvector = [ 0.73725490 0.82352941 0.93333333 ];
	case 'LightSteelBlue3',
	 colorvector = [ 0.63529412 0.70980392 0.80392157 ];
	case 'LightSteelBlue4',
	 colorvector = [ 0.43137255 0.48235294 0.54509804 ];
	case 'LightBlue1',
	 colorvector = [ 0.74901961 0.93725490 1 ];
	case 'LightBlue2',
	 colorvector = [ 0.69803922 0.87450980 0.93333333 ];
	case 'LightBlue3',
	 colorvector = [ 0.60392157 0.75294118 0.80392157 ];
	case 'LightBlue4',
	 colorvector = [ 0.40784314 0.51372549 0.54509804 ];
	case 'LightCyan1',
	 colorvector = [ 0.87843137 1 1 ];
	case 'LightCyan2',
	 colorvector = [ 0.81960784 0.93333333 0.93333333 ];
	case 'LightCyan3',
	 colorvector = [ 0.70588235 0.80392157 0.80392157 ];
	case 'LightCyan4',
	 colorvector = [ 0.47843137 0.54509804 0.54509804 ];
	case 'PaleTurquoise1',
	 colorvector = [ 0.73333333 1 1 ];
	case 'PaleTurquoise2',
	 colorvector = [ 0.68235294 0.93333333 0.93333333 ];
	case 'PaleTurquoise3',
	 colorvector = [ 0.58823529 0.80392157 0.80392157 ];
	case 'PaleTurquoise4',
	 colorvector = [ 0.4 0.54509804 0.54509804 ];
	case 'CadetBlue1',
	 colorvector = [ 0.59607843 0.96078431 1 ];
	case 'CadetBlue2',
	 colorvector = [ 0.55686275 0.89803922 0.93333333 ];
	case 'CadetBlue3',
	 colorvector = [ 0.47843137 0.77254902 0.80392157 ];
	case 'CadetBlue4',
	 colorvector = [ 0.32549020 0.52549020 0.54509804 ];
	case 'turquoise1',
	 colorvector = [ 0 0.96078431 1 ];
	case 'turquoise2',
	 colorvector = [ 0 0.89803922 0.93333333 ];
	case 'turquoise3',
	 colorvector = [ 0 0.77254902 0.80392157 ];
	case 'turquoise4',
	 colorvector = [ 0 0.52549020 0.54509804 ];
	case 'cyan1',
	 colorvector = [ 0 1 1 ];
	case 'cyan2',
	 colorvector = [ 0 0.93333333 0.93333333 ];
	case 'cyan3',
	 colorvector = [ 0 0.80392157 0.80392157 ];
	case 'cyan4',
	 colorvector = [ 0 0.54509804 0.54509804 ];
	case 'DarkSlateGray1',
	 colorvector = [ 0.59215686 1 1 ];
	case 'DarkSlateGray2',
	 colorvector = [ 0.55294118 0.93333333 0.93333333 ];
	case 'DarkSlateGray3',
	 colorvector = [ 0.47450980 0.80392157 0.80392157 ];
	case 'DarkSlateGray4',
	 colorvector = [ 0.32156863 0.54509804 0.54509804 ];
	case 'aquamarine1',
	 colorvector = [ 0.49803922 1 0.83137255 ];
	case 'aquamarine2',
	 colorvector = [ 0.46274510 0.93333333 0.77647059 ];
	case 'aquamarine3',
	 colorvector = [ 0.4 0.80392157 0.66666667 ];
	case 'aquamarine4',
	 colorvector = [ 0.27058824 0.54509804 0.45490196 ];
	case 'DarkSeaGreen1',
	 colorvector = [ 0.75686275 1 0.75686275 ];
	case 'DarkSeaGreen2',
	 colorvector = [ 0.70588235 0.93333333 0.70588235 ];
	case 'DarkSeaGreen3',
	 colorvector = [ 0.60784314 0.80392157 0.60784314 ];
	case 'DarkSeaGreen4',
	 colorvector = [ 0.41176471 0.54509804 0.41176471 ];
	case 'SeaGreen1',
	 colorvector = [ 0.32941176 1 0.62352941 ];
	case 'SeaGreen2',
	 colorvector = [ 0.30588235 0.93333333 0.58039216 ];
	case 'SeaGreen3',
	 colorvector = [ 0.26274510 0.80392157 0.50196078 ];
	case 'SeaGreen4',
	 colorvector = [ 0.18039216 0.54509804 0.34117647 ];
	case 'PaleGreen1',
	 colorvector = [ 0.60392157 1 0.60392157 ];
	case 'PaleGreen2',
	 colorvector = [ 0.56470588 0.93333333 0.56470588 ];
	case 'PaleGreen3',
	 colorvector = [ 0.48627451 0.80392157 0.48627451 ];
	case 'PaleGreen4',
	 colorvector = [ 0.32941176 0.54509804 0.32941176 ];
	case 'SpringGreen1',
	 colorvector = [ 0 1 0.49803922 ];
	case 'SpringGreen2',
	 colorvector = [ 0 0.93333333 0.46274510 ];
	case 'SpringGreen3',
	 colorvector = [ 0 0.80392157 0.4 ];
	case 'SpringGreen4',
	 colorvector = [ 0 0.54509804 0.27058824 ];
	case 'green1',
	 colorvector = [ 0 1 0 ];
	case 'green2',
	 colorvector = [ 0 0.93333333 0 ];
	case 'green3',
	 colorvector = [ 0 0.80392157 0 ];
	case 'green4',
	 colorvector = [ 0 0.54509804 0 ];
	case 'chartreuse1',
	 colorvector = [ 0.49803922 1 0 ];
	case 'chartreuse2',
	 colorvector = [ 0.46274510 0.93333333 0 ];
	case 'chartreuse3',
	 colorvector = [ 0.4 0.80392157 0 ];
	case 'chartreuse4',
	 colorvector = [ 0.27058824 0.54509804 0 ];
	case 'OliveDrab1',
	 colorvector = [ 0.75294118 1 0.24313725 ];
	case 'OliveDrab2',
	 colorvector = [ 0.70196078 0.93333333 0.22745098 ];
	case 'OliveDrab3',
	 colorvector = [ 0.60392157 0.80392157 0.19607843 ];
	case 'OliveDrab4',
	 colorvector = [ 0.41176471 0.54509804 0.13333333 ];
	case 'DarkOliveGreen1',
	 colorvector = [ 0.79215686 1 0.43921569 ];
	case 'DarkOliveGreen2',
	 colorvector = [ 0.73725490 0.93333333 0.40784314 ];
	case 'DarkOliveGreen3',
	 colorvector = [ 0.63529412 0.80392157 0.35294118 ];
	case 'DarkOliveGreen4',
	 colorvector = [ 0.43137255 0.54509804 0.23921569 ];
	case 'khaki1',
	 colorvector = [ 1 0.96470588 0.56078431 ];
	case 'khaki2',
	 colorvector = [ 0.93333333 0.90196078 0.52156863 ];
	case 'khaki3',
	 colorvector = [ 0.80392157 0.77647059 0.45098039 ];
	case 'khaki4',
	 colorvector = [ 0.54509804 0.52549020 0.30588235 ];
	case 'LightGoldenrod1',
	 colorvector = [ 1 0.92549020 0.54509804 ];
	case 'LightGoldenrod2',
	 colorvector = [ 0.93333333 0.86274510 0.50980392 ];
	case 'LightGoldenrod3',
	 colorvector = [ 0.80392157 0.74509804 0.43921569 ];
	case 'LightGoldenrod4',
	 colorvector = [ 0.54509804 0.50588235 0.29803922 ];
	case 'LightYellow1',
	 colorvector = [ 1 1 0.87843137 ];
	case 'LightYellow2',
	 colorvector = [ 0.93333333 0.93333333 0.81960784 ];
	case 'LightYellow3',
	 colorvector = [ 0.80392157 0.80392157 0.70588235 ];
	case 'LightYellow4',
	 colorvector = [ 0.54509804 0.54509804 0.47843137 ];
	case 'yellow1',
	 colorvector = [ 1 1 0 ];
	case 'yellow2',
	 colorvector = [ 0.93333333 0.93333333 0 ];
	case 'yellow3',
	 colorvector = [ 0.80392157 0.80392157 0 ];
	case 'yellow4',
	 colorvector = [ 0.54509804 0.54509804 0 ];
	case 'gold1',
	 colorvector = [ 1 0.84313725 0 ];
	case 'gold2',
	 colorvector = [ 0.93333333 0.78823529 0 ];
	case 'gold3',
	 colorvector = [ 0.80392157 0.67843137 0 ];
	case 'gold4',
	 colorvector = [ 0.54509804 0.45882353 0 ];
	case 'goldenrod1',
	 colorvector = [ 1 0.75686275 0.14509804 ];
	case 'goldenrod2',
	 colorvector = [ 0.93333333 0.70588235 0.13333333 ];
	case 'goldenrod3',
	 colorvector = [ 0.80392157 0.60784314 0.11372549 ];
	case 'goldenrod4',
	 colorvector = [ 0.54509804 0.41176471 0.078431373 ];
	case 'DarkGoldenrod1',
	 colorvector = [ 1 0.72549020 0.058823529 ];
	case 'DarkGoldenrod2',
	 colorvector = [ 0.93333333 0.67843137 0.054901961 ];
	case 'DarkGoldenrod3',
	 colorvector = [ 0.80392157 0.58431373 0.047058824 ];
	case 'DarkGoldenrod4',
	 colorvector = [ 0.54509804 0.39607843 0.031372549 ];
	case 'RosyBrown1',
	 colorvector = [ 1 0.75686275 0.75686275 ];
	case 'RosyBrown2',
	 colorvector = [ 0.93333333 0.70588235 0.70588235 ];
	case 'RosyBrown3',
	 colorvector = [ 0.80392157 0.60784314 0.60784314 ];
	case 'RosyBrown4',
	 colorvector = [ 0.54509804 0.41176471 0.41176471 ];
	case 'IndianRed1',
	 colorvector = [ 1 0.41568627 0.41568627 ];
	case 'IndianRed2',
	 colorvector = [ 0.93333333 0.38823529 0.38823529 ];
	case 'IndianRed3',
	 colorvector = [ 0.80392157 0.33333333 0.33333333 ];
	case 'IndianRed4',
	 colorvector = [ 0.54509804 0.22745098 0.22745098 ];
	case 'sienna1',
	 colorvector = [ 1 0.50980392 0.27843137 ];
	case 'sienna2',
	 colorvector = [ 0.93333333 0.47450980 0.25882353 ];
	case 'sienna3',
	 colorvector = [ 0.80392157 0.40784314 0.22352941 ];
	case 'sienna4',
	 colorvector = [ 0.54509804 0.27843137 0.14901961 ];
	case 'burlywood1',
	 colorvector = [ 1 0.82745098 0.60784314 ];
	case 'burlywood2',
	 colorvector = [ 0.93333333 0.77254902 0.56862745 ];
	case 'burlywood3',
	 colorvector = [ 0.80392157 0.66666667 0.49019608 ];
	case 'burlywood4',
	 colorvector = [ 0.54509804 0.45098039 0.33333333 ];
	case 'wheat1',
	 colorvector = [ 1 0.90588235 0.72941176 ];
	case 'wheat2',
	 colorvector = [ 0.93333333 0.84705882 0.68235294 ];
	case 'wheat3',
	 colorvector = [ 0.80392157 0.72941176 0.58823529 ];
	case 'wheat4',
	 colorvector = [ 0.54509804 0.49411765 0.4 ];
	case 'tan1',
	 colorvector = [ 1 0.64705882 0.30980392 ];
	case 'tan2',
	 colorvector = [ 0.93333333 0.60392157 0.28627451 ];
	case 'tan3',
	 colorvector = [ 0.80392157 0.52156863 0.24705882 ];
	case 'tan4',
	 colorvector = [ 0.54509804 0.35294118 0.16862745 ];
	case 'chocolate1',
	 colorvector = [ 1 0.49803922 0.14117647 ];
	case 'chocolate2',
	 colorvector = [ 0.93333333 0.46274510 0.12941176 ];
	case 'chocolate3',
	 colorvector = [ 0.80392157 0.4 0.11372549 ];
	case 'chocolate4',
	 colorvector = [ 0.54509804 0.27058824 0.074509804 ];
	case 'firebrick1',
	 colorvector = [ 1 0.18823529 0.18823529 ];
	case 'firebrick2',
	 colorvector = [ 0.93333333 0.17254902 0.17254902 ];
	case 'firebrick3',
	 colorvector = [ 0.80392157 0.14901961 0.14901961 ];
	case 'firebrick4',
	 colorvector = [ 0.54509804 0.10196078 0.10196078 ];
	case 'brown1',
	 colorvector = [ 1 0.25098039 0.25098039 ];
	case 'brown2',
	 colorvector = [ 0.93333333 0.23137255 0.23137255 ];
	case 'brown3',
	 colorvector = [ 0.80392157 0.2 0.2 ];
	case 'brown4',
	 colorvector = [ 0.54509804 0.13725490 0.13725490 ];
	case 'salmon1',
	 colorvector = [ 1 0.54901961 0.41176471 ];
	case 'salmon2',
	 colorvector = [ 0.93333333 0.50980392 0.38431373 ];
	case 'salmon3',
	 colorvector = [ 0.80392157 0.43921569 0.32941176 ];
	case 'salmon4',
	 colorvector = [ 0.54509804 0.29803922 0.22352941 ];
	case 'LightSalmon1',
	 colorvector = [ 1 0.62745098 0.47843137 ];
	case 'LightSalmon2',
	 colorvector = [ 0.93333333 0.58431373 0.44705882 ];
	case 'LightSalmon3',
	 colorvector = [ 0.80392157 0.50588235 0.38431373 ];
	case 'LightSalmon4',
	 colorvector = [ 0.54509804 0.34117647 0.25882353 ];
	case 'orange1',
	 colorvector = [ 1 0.64705882 0 ];
	case 'orange2',
	 colorvector = [ 0.93333333 0.60392157 0 ];
	case 'orange3',
	 colorvector = [ 0.80392157 0.52156863 0 ];
	case 'orange4',
	 colorvector = [ 0.54509804 0.35294118 0 ];
	case 'DarkOrange1',
	 colorvector = [ 1 0.49803922 0 ];
	case 'DarkOrange2',
	 colorvector = [ 0.93333333 0.46274510 0 ];
	case 'DarkOrange3',
	 colorvector = [ 0.80392157 0.4 0 ];
	case 'DarkOrange4',
	 colorvector = [ 0.54509804 0.27058824 0 ];
	case 'coral1',
	 colorvector = [ 1 0.44705882 0.33725490 ];
	case 'coral2',
	 colorvector = [ 0.93333333 0.41568627 0.31372549 ];
	case 'coral3',
	 colorvector = [ 0.80392157 0.35686275 0.27058824 ];
	case 'coral4',
	 colorvector = [ 0.54509804 0.24313725 0.18431373 ];
	case 'tomato1',
	 colorvector = [ 1 0.38823529 0.27843137 ];
	case 'tomato2',
	 colorvector = [ 0.93333333 0.36078431 0.25882353 ];
	case 'tomato3',
	 colorvector = [ 0.80392157 0.30980392 0.22352941 ];
	case 'tomato4',
	 colorvector = [ 0.54509804 0.21176471 0.14901961 ];
	case 'OrangeRed1',
	 colorvector = [ 1 0.27058824 0 ];
	case 'OrangeRed2',
	 colorvector = [ 0.93333333 0.25098039 0 ];
	case 'OrangeRed3',
	 colorvector = [ 0.80392157 0.21568627 0 ];
	case 'OrangeRed4',
	 colorvector = [ 0.54509804 0.14509804 0 ];
	case 'red1',
	 colorvector = [ 1 0 0 ];
	case 'red2',
	 colorvector = [ 0.93333333 0 0 ];
	case 'red3',
	 colorvector = [ 0.80392157 0 0 ];
	case 'red4',
	 colorvector = [ 0.54509804 0 0 ];
	case 'DeepPink1',
	 colorvector = [ 1 0.078431373 0.57647059 ];
	case 'DeepPink2',
	 colorvector = [ 0.93333333 0.070588235 0.53725490 ];
	case 'DeepPink3',
	 colorvector = [ 0.80392157 0.062745098 0.46274510 ];
	case 'DeepPink4',
	 colorvector = [ 0.54509804 0.039215686 0.31372549 ];
	case 'HotPink1',
	 colorvector = [ 1 0.43137255 0.70588235 ];
	case 'HotPink2',
	 colorvector = [ 0.93333333 0.41568627 0.65490196 ];
	case 'HotPink3',
	 colorvector = [ 0.80392157 0.37647059 0.56470588 ];
	case 'HotPink4',
	 colorvector = [ 0.54509804 0.22745098 0.38431373 ];
	case 'pink1',
	 colorvector = [ 1 0.70980392 0.77254902 ];
	case 'pink2',
	 colorvector = [ 0.93333333 0.66274510 0.72156863 ];
	case 'pink3',
	 colorvector = [ 0.80392157 0.56862745 0.61960784 ];
	case 'pink4',
	 colorvector = [ 0.54509804 0.38823529 0.42352941 ];
	case 'LightPink1',
	 colorvector = [ 1 0.68235294 0.72549020 ];
	case 'LightPink2',
	 colorvector = [ 0.93333333 0.63529412 0.67843137 ];
	case 'LightPink3',
	 colorvector = [ 0.80392157 0.54901961 0.58431373 ];
	case 'LightPink4',
	 colorvector = [ 0.54509804 0.37254902 0.39607843 ];
	case 'PaleVioletRed1',
	 colorvector = [ 1 0.50980392 0.67058824 ];
	case 'PaleVioletRed2',
	 colorvector = [ 0.93333333 0.47450980 0.62352941 ];
	case 'PaleVioletRed3',
	 colorvector = [ 0.80392157 0.40784314 0.53725490 ];
	case 'PaleVioletRed4',
	 colorvector = [ 0.54509804 0.27843137 0.36470588 ];
	case 'maroon1',
	 colorvector = [ 1 0.20392157 0.70196078 ];
	case 'maroon2',
	 colorvector = [ 0.93333333 0.18823529 0.65490196 ];
	case 'maroon3',
	 colorvector = [ 0.80392157 0.16078431 0.56470588 ];
	case 'maroon4',
	 colorvector = [ 0.54509804 0.10980392 0.38431373 ];
	case 'VioletRed1',
	 colorvector = [ 1 0.24313725 0.58823529 ];
	case 'VioletRed2',
	 colorvector = [ 0.93333333 0.22745098 0.54901961 ];
	case 'VioletRed3',
	 colorvector = [ 0.80392157 0.19607843 0.47058824 ];
	case 'VioletRed4',
	 colorvector = [ 0.54509804 0.13333333 0.32156863 ];
	case 'magenta1',
	 colorvector = [ 1 0 1 ];
	case 'magenta2',
	 colorvector = [ 0.93333333 0 0.93333333 ];
	case 'magenta3',
	 colorvector = [ 0.80392157 0 0.80392157 ];
	case 'magenta4',
	 colorvector = [ 0.54509804 0 0.54509804 ];
	case 'orchid1',
	 colorvector = [ 1 0.51372549 0.98039216 ];
	case 'orchid2',
	 colorvector = [ 0.93333333 0.47843137 0.91372549 ];
	case 'orchid3',
	 colorvector = [ 0.80392157 0.41176471 0.78823529 ];
	case 'orchid4',
	 colorvector = [ 0.54509804 0.27843137 0.53725490 ];
	case 'plum1',
	 colorvector = [ 1 0.73333333 1 ];
	case 'plum2',
	 colorvector = [ 0.93333333 0.68235294 0.93333333 ];
	case 'plum3',
	 colorvector = [ 0.80392157 0.58823529 0.80392157 ];
	case 'plum4',
	 colorvector = [ 0.54509804 0.4 0.54509804 ];
	case 'MediumOrchid1',
	 colorvector = [ 0.87843137 0.4 1 ];
	case 'MediumOrchid2',
	 colorvector = [ 0.81960784 0.37254902 0.93333333 ];
	case 'MediumOrchid3',
	 colorvector = [ 0.70588235 0.32156863 0.80392157 ];
	case 'MediumOrchid4',
	 colorvector = [ 0.47843137 0.21568627 0.54509804 ];
	case 'DarkOrchid1',
	 colorvector = [ 0.74901961 0.24313725 1 ];
	case 'DarkOrchid2',
	 colorvector = [ 0.69803922 0.22745098 0.93333333 ];
	case 'DarkOrchid3',
	 colorvector = [ 0.60392157 0.19607843 0.80392157 ];
	case 'DarkOrchid4',
	 colorvector = [ 0.40784314 0.13333333 0.54509804 ];
	case 'purple1',
	 colorvector = [ 0.60784314 0.18823529 1 ];
	case 'purple2',
	 colorvector = [ 0.56862745 0.17254902 0.93333333 ];
	case 'purple3',
	 colorvector = [ 0.49019608 0.14901961 0.80392157 ];
	case 'purple4',
	 colorvector = [ 0.33333333 0.10196078 0.54509804 ];
	case 'MediumPurple1',
	 colorvector = [ 0.67058824 0.50980392 1 ];
	case 'MediumPurple2',
	 colorvector = [ 0.62352941 0.47450980 0.93333333 ];
	case 'MediumPurple3',
	 colorvector = [ 0.53725490 0.40784314 0.80392157 ];
	case 'MediumPurple4',
	 colorvector = [ 0.36470588 0.27843137 0.54509804 ];
	case 'thistle1',
	 colorvector = [ 1 0.88235294 1 ];
	case 'thistle2',
	 colorvector = [ 0.93333333 0.82352941 0.93333333 ];
	case 'thistle3',
	 colorvector = [ 0.80392157 0.70980392 0.80392157 ];
	case 'thistle4',
	 colorvector = [ 0.54509804 0.48235294 0.54509804 ];
	case 'gray0',
	 colorvector = [ 0 0 0 ];
	case 'grey0',
	 colorvector = [ 0 0 0 ];
	case 'gray1',
	 colorvector = [ 0.011764706 0.011764706 0.011764706 ];
	case 'grey1',
	 colorvector = [ 0.011764706 0.011764706 0.011764706 ];
	case 'gray2',
	 colorvector = [ 0.019607843 0.019607843 0.019607843 ];
	case 'grey2',
	 colorvector = [ 0.019607843 0.019607843 0.019607843 ];
	case 'gray3',
	 colorvector = [ 0.031372549 0.031372549 0.031372549 ];
	case 'grey3',
	 colorvector = [ 0.031372549 0.031372549 0.031372549 ];
	case 'gray4',
	 colorvector = [ 0.039215686 0.039215686 0.039215686 ];
	case 'grey4',
	 colorvector = [ 0.039215686 0.039215686 0.039215686 ];
	case 'gray5',
	 colorvector = [ 0.050980392 0.050980392 0.050980392 ];
	case 'grey5',
	 colorvector = [ 0.050980392 0.050980392 0.050980392 ];
	case 'gray6',
	 colorvector = [ 0.058823529 0.058823529 0.058823529 ];
	case 'grey6',
	 colorvector = [ 0.058823529 0.058823529 0.058823529 ];
	case 'gray7',
	 colorvector = [ 0.070588235 0.070588235 0.070588235 ];
	case 'grey7',
	 colorvector = [ 0.070588235 0.070588235 0.070588235 ];
	case 'gray8',
	 colorvector = [ 0.078431373 0.078431373 0.078431373 ];
	case 'grey8',
	 colorvector = [ 0.078431373 0.078431373 0.078431373 ];
	case 'gray9',
	 colorvector = [ 0.090196078 0.090196078 0.090196078 ];
	case 'grey9',
	 colorvector = [ 0.090196078 0.090196078 0.090196078 ];
	case 'gray10',
	 colorvector = [ 0.10196078 0.10196078 0.10196078 ];
	case 'grey10',
	 colorvector = [ 0.10196078 0.10196078 0.10196078 ];
	case 'gray11',
	 colorvector = [ 0.10980392 0.10980392 0.10980392 ];
	case 'grey11',
	 colorvector = [ 0.10980392 0.10980392 0.10980392 ];
	case 'gray12',
	 colorvector = [ 0.12156863 0.12156863 0.12156863 ];
	case 'grey12',
	 colorvector = [ 0.12156863 0.12156863 0.12156863 ];
	case 'gray13',
	 colorvector = [ 0.12941176 0.12941176 0.12941176 ];
	case 'grey13',
	 colorvector = [ 0.12941176 0.12941176 0.12941176 ];
	case 'gray14',
	 colorvector = [ 0.14117647 0.14117647 0.14117647 ];
	case 'grey14',
	 colorvector = [ 0.14117647 0.14117647 0.14117647 ];
	case 'gray15',
	 colorvector = [ 0.14901961 0.14901961 0.14901961 ];
	case 'grey15',
	 colorvector = [ 0.14901961 0.14901961 0.14901961 ];
	case 'gray16',
	 colorvector = [ 0.16078431 0.16078431 0.16078431 ];
	case 'grey16',
	 colorvector = [ 0.16078431 0.16078431 0.16078431 ];
	case 'gray17',
	 colorvector = [ 0.16862745 0.16862745 0.16862745 ];
	case 'grey17',
	 colorvector = [ 0.16862745 0.16862745 0.16862745 ];
	case 'gray18',
	 colorvector = [ 0.18039216 0.18039216 0.18039216 ];
	case 'grey18',
	 colorvector = [ 0.18039216 0.18039216 0.18039216 ];
	case 'gray19',
	 colorvector = [ 0.18823529 0.18823529 0.18823529 ];
	case 'grey19',
	 colorvector = [ 0.18823529 0.18823529 0.18823529 ];
	case 'gray20',
	 colorvector = [ 0.2 0.2 0.2 ];
	case 'grey20',
	 colorvector = [ 0.2 0.2 0.2 ];
	case 'gray21',
	 colorvector = [ 0.21176471 0.21176471 0.21176471 ];
	case 'grey21',
	 colorvector = [ 0.21176471 0.21176471 0.21176471 ];
	case 'gray22',
	 colorvector = [ 0.21960784 0.21960784 0.21960784 ];
	case 'grey22',
	 colorvector = [ 0.21960784 0.21960784 0.21960784 ];
	case 'gray23',
	 colorvector = [ 0.23137255 0.23137255 0.23137255 ];
	case 'grey23',
	 colorvector = [ 0.23137255 0.23137255 0.23137255 ];
	case 'gray24',
	 colorvector = [ 0.23921569 0.23921569 0.23921569 ];
	case 'grey24',
	 colorvector = [ 0.23921569 0.23921569 0.23921569 ];
	case 'gray25',
	 colorvector = [ 0.25098039 0.25098039 0.25098039 ];
	case 'grey25',
	 colorvector = [ 0.25098039 0.25098039 0.25098039 ];
	case 'gray26',
	 colorvector = [ 0.25882353 0.25882353 0.25882353 ];
	case 'grey26',
	 colorvector = [ 0.25882353 0.25882353 0.25882353 ];
	case 'gray27',
	 colorvector = [ 0.27058824 0.27058824 0.27058824 ];
	case 'grey27',
	 colorvector = [ 0.27058824 0.27058824 0.27058824 ];
	case 'gray28',
	 colorvector = [ 0.27843137 0.27843137 0.27843137 ];
	case 'grey28',
	 colorvector = [ 0.27843137 0.27843137 0.27843137 ];
	case 'gray29',
	 colorvector = [ 0.29019608 0.29019608 0.29019608 ];
	case 'grey29',
	 colorvector = [ 0.29019608 0.29019608 0.29019608 ];
	case 'gray30',
	 colorvector = [ 0.30196078 0.30196078 0.30196078 ];
	case 'grey30',
	 colorvector = [ 0.30196078 0.30196078 0.30196078 ];
	case 'gray31',
	 colorvector = [ 0.30980392 0.30980392 0.30980392 ];
	case 'grey31',
	 colorvector = [ 0.30980392 0.30980392 0.30980392 ];
	case 'gray32',
	 colorvector = [ 0.32156863 0.32156863 0.32156863 ];
	case 'grey32',
	 colorvector = [ 0.32156863 0.32156863 0.32156863 ];
	case 'gray33',
	 colorvector = [ 0.32941176 0.32941176 0.32941176 ];
	case 'grey33',
	 colorvector = [ 0.32941176 0.32941176 0.32941176 ];
	case 'gray34',
	 colorvector = [ 0.34117647 0.34117647 0.34117647 ];
	case 'grey34',
	 colorvector = [ 0.34117647 0.34117647 0.34117647 ];
	case 'gray35',
	 colorvector = [ 0.34901961 0.34901961 0.34901961 ];
	case 'grey35',
	 colorvector = [ 0.34901961 0.34901961 0.34901961 ];
	case 'gray36',
	 colorvector = [ 0.36078431 0.36078431 0.36078431 ];
	case 'grey36',
	 colorvector = [ 0.36078431 0.36078431 0.36078431 ];
	case 'gray37',
	 colorvector = [ 0.36862745 0.36862745 0.36862745 ];
	case 'grey37',
	 colorvector = [ 0.36862745 0.36862745 0.36862745 ];
	case 'gray38',
	 colorvector = [ 0.38039216 0.38039216 0.38039216 ];
	case 'grey38',
	 colorvector = [ 0.38039216 0.38039216 0.38039216 ];
	case 'gray39',
	 colorvector = [ 0.38823529 0.38823529 0.38823529 ];
	case 'grey39',
	 colorvector = [ 0.38823529 0.38823529 0.38823529 ];
	case 'gray40',
	 colorvector = [ 0.4 0.4 0.4 ];
	case 'grey40',
	 colorvector = [ 0.4 0.4 0.4 ];
	case 'gray41',
	 colorvector = [ 0.41176471 0.41176471 0.41176471 ];
	case 'grey41',
	 colorvector = [ 0.41176471 0.41176471 0.41176471 ];
	case 'gray42',
	 colorvector = [ 0.41960784 0.41960784 0.41960784 ];
	case 'grey42',
	 colorvector = [ 0.41960784 0.41960784 0.41960784 ];
	case 'gray43',
	 colorvector = [ 0.43137255 0.43137255 0.43137255 ];
	case 'grey43',
	 colorvector = [ 0.43137255 0.43137255 0.43137255 ];
	case 'gray44',
	 colorvector = [ 0.43921569 0.43921569 0.43921569 ];
	case 'grey44',
	 colorvector = [ 0.43921569 0.43921569 0.43921569 ];
	case 'gray45',
	 colorvector = [ 0.45098039 0.45098039 0.45098039 ];
	case 'grey45',
	 colorvector = [ 0.45098039 0.45098039 0.45098039 ];
	case 'gray46',
	 colorvector = [ 0.45882353 0.45882353 0.45882353 ];
	case 'grey46',
	 colorvector = [ 0.45882353 0.45882353 0.45882353 ];
	case 'gray47',
	 colorvector = [ 0.47058824 0.47058824 0.47058824 ];
	case 'grey47',
	 colorvector = [ 0.47058824 0.47058824 0.47058824 ];
	case 'gray48',
	 colorvector = [ 0.47843137 0.47843137 0.47843137 ];
	case 'grey48',
	 colorvector = [ 0.47843137 0.47843137 0.47843137 ];
	case 'gray49',
	 colorvector = [ 0.49019608 0.49019608 0.49019608 ];
	case 'grey49',
	 colorvector = [ 0.49019608 0.49019608 0.49019608 ];
	case 'gray50',
	 colorvector = [ 0.49803922 0.49803922 0.49803922 ];
	case 'grey50',
	 colorvector = [ 0.49803922 0.49803922 0.49803922 ];
	case 'gray51',
	 colorvector = [ 0.50980392 0.50980392 0.50980392 ];
	case 'grey51',
	 colorvector = [ 0.50980392 0.50980392 0.50980392 ];
	case 'gray52',
	 colorvector = [ 0.52156863 0.52156863 0.52156863 ];
	case 'grey52',
	 colorvector = [ 0.52156863 0.52156863 0.52156863 ];
	case 'gray53',
	 colorvector = [ 0.52941176 0.52941176 0.52941176 ];
	case 'grey53',
	 colorvector = [ 0.52941176 0.52941176 0.52941176 ];
	case 'gray54',
	 colorvector = [ 0.54117647 0.54117647 0.54117647 ];
	case 'grey54',
	 colorvector = [ 0.54117647 0.54117647 0.54117647 ];
	case 'gray55',
	 colorvector = [ 0.54901961 0.54901961 0.54901961 ];
	case 'grey55',
	 colorvector = [ 0.54901961 0.54901961 0.54901961 ];
	case 'gray56',
	 colorvector = [ 0.56078431 0.56078431 0.56078431 ];
	case 'grey56',
	 colorvector = [ 0.56078431 0.56078431 0.56078431 ];
	case 'gray57',
	 colorvector = [ 0.56862745 0.56862745 0.56862745 ];
	case 'grey57',
	 colorvector = [ 0.56862745 0.56862745 0.56862745 ];
	case 'gray58',
	 colorvector = [ 0.58039216 0.58039216 0.58039216 ];
	case 'grey58',
	 colorvector = [ 0.58039216 0.58039216 0.58039216 ];
	case 'gray59',
	 colorvector = [ 0.58823529 0.58823529 0.58823529 ];
	case 'grey59',
	 colorvector = [ 0.58823529 0.58823529 0.58823529 ];
	case 'gray60',
	 colorvector = [ 0.6 0.6 0.6 ];
	case 'grey60',
	 colorvector = [ 0.6 0.6 0.6 ];
	case 'gray61',
	 colorvector = [ 0.61176471 0.61176471 0.61176471 ];
	case 'grey61',
	 colorvector = [ 0.61176471 0.61176471 0.61176471 ];
	case 'gray62',
	 colorvector = [ 0.61960784 0.61960784 0.61960784 ];
	case 'grey62',
	 colorvector = [ 0.61960784 0.61960784 0.61960784 ];
	case 'gray63',
	 colorvector = [ 0.63137255 0.63137255 0.63137255 ];
	case 'grey63',
	 colorvector = [ 0.63137255 0.63137255 0.63137255 ];
	case 'gray64',
	 colorvector = [ 0.63921569 0.63921569 0.63921569 ];
	case 'grey64',
	 colorvector = [ 0.63921569 0.63921569 0.63921569 ];
	case 'gray65',
	 colorvector = [ 0.65098039 0.65098039 0.65098039 ];
	case 'grey65',
	 colorvector = [ 0.65098039 0.65098039 0.65098039 ];
	case 'gray66',
	 colorvector = [ 0.65882353 0.65882353 0.65882353 ];
	case 'grey66',
	 colorvector = [ 0.65882353 0.65882353 0.65882353 ];
	case 'gray67',
	 colorvector = [ 0.67058824 0.67058824 0.67058824 ];
	case 'grey67',
	 colorvector = [ 0.67058824 0.67058824 0.67058824 ];
	case 'gray68',
	 colorvector = [ 0.67843137 0.67843137 0.67843137 ];
	case 'grey68',
	 colorvector = [ 0.67843137 0.67843137 0.67843137 ];
	case 'gray69',
	 colorvector = [ 0.69019608 0.69019608 0.69019608 ];
	case 'grey69',
	 colorvector = [ 0.69019608 0.69019608 0.69019608 ];
	case 'gray70',
	 colorvector = [ 0.70196078 0.70196078 0.70196078 ];
	case 'grey70',
	 colorvector = [ 0.70196078 0.70196078 0.70196078 ];
	case 'gray71',
	 colorvector = [ 0.70980392 0.70980392 0.70980392 ];
	case 'grey71',
	 colorvector = [ 0.70980392 0.70980392 0.70980392 ];
	case 'gray72',
	 colorvector = [ 0.72156863 0.72156863 0.72156863 ];
	case 'grey72',
	 colorvector = [ 0.72156863 0.72156863 0.72156863 ];
	case 'gray73',
	 colorvector = [ 0.72941176 0.72941176 0.72941176 ];
	case 'grey73',
	 colorvector = [ 0.72941176 0.72941176 0.72941176 ];
	case 'gray74',
	 colorvector = [ 0.74117647 0.74117647 0.74117647 ];
	case 'grey74',
	 colorvector = [ 0.74117647 0.74117647 0.74117647 ];
	case 'gray75',
	 colorvector = [ 0.74901961 0.74901961 0.74901961 ];
	case 'grey75',
	 colorvector = [ 0.74901961 0.74901961 0.74901961 ];
	case 'gray76',
	 colorvector = [ 0.76078431 0.76078431 0.76078431 ];
	case 'grey76',
	 colorvector = [ 0.76078431 0.76078431 0.76078431 ];
	case 'gray77',
	 colorvector = [ 0.76862745 0.76862745 0.76862745 ];
	case 'grey77',
	 colorvector = [ 0.76862745 0.76862745 0.76862745 ];
	case 'gray78',
	 colorvector = [ 0.78039216 0.78039216 0.78039216 ];
	case 'grey78',
	 colorvector = [ 0.78039216 0.78039216 0.78039216 ];
	case 'gray79',
	 colorvector = [ 0.78823529 0.78823529 0.78823529 ];
	case 'grey79',
	 colorvector = [ 0.78823529 0.78823529 0.78823529 ];
	case 'gray80',
	 colorvector = [ 0.8 0.8 0.8 ];
	case 'grey80',
	 colorvector = [ 0.8 0.8 0.8 ];
	case 'gray81',
	 colorvector = [ 0.81176471 0.81176471 0.81176471 ];
	case 'grey81',
	 colorvector = [ 0.81176471 0.81176471 0.81176471 ];
	case 'gray82',
	 colorvector = [ 0.81960784 0.81960784 0.81960784 ];
	case 'grey82',
	 colorvector = [ 0.81960784 0.81960784 0.81960784 ];
	case 'gray83',
	 colorvector = [ 0.83137255 0.83137255 0.83137255 ];
	case 'grey83',
	 colorvector = [ 0.83137255 0.83137255 0.83137255 ];
	case 'gray84',
	 colorvector = [ 0.83921569 0.83921569 0.83921569 ];
	case 'grey84',
	 colorvector = [ 0.83921569 0.83921569 0.83921569 ];
	case 'gray85',
	 colorvector = [ 0.85098039 0.85098039 0.85098039 ];
	case 'grey85',
	 colorvector = [ 0.85098039 0.85098039 0.85098039 ];
	case 'gray86',
	 colorvector = [ 0.85882353 0.85882353 0.85882353 ];
	case 'grey86',
	 colorvector = [ 0.85882353 0.85882353 0.85882353 ];
	case 'gray87',
	 colorvector = [ 0.87058824 0.87058824 0.87058824 ];
	case 'grey87',
	 colorvector = [ 0.87058824 0.87058824 0.87058824 ];
	case 'gray88',
	 colorvector = [ 0.87843137 0.87843137 0.87843137 ];
	case 'grey88',
	 colorvector = [ 0.87843137 0.87843137 0.87843137 ];
	case 'gray89',
	 colorvector = [ 0.89019608 0.89019608 0.89019608 ];
	case 'grey89',
	 colorvector = [ 0.89019608 0.89019608 0.89019608 ];
	case 'gray90',
	 colorvector = [ 0.89803922 0.89803922 0.89803922 ];
	case 'grey90',
	 colorvector = [ 0.89803922 0.89803922 0.89803922 ];
	case 'gray91',
	 colorvector = [ 0.90980392 0.90980392 0.90980392 ];
	case 'grey91',
	 colorvector = [ 0.90980392 0.90980392 0.90980392 ];
	case 'gray92',
	 colorvector = [ 0.92156863 0.92156863 0.92156863 ];
	case 'grey92',
	 colorvector = [ 0.92156863 0.92156863 0.92156863 ];
	case 'gray93',
	 colorvector = [ 0.92941176 0.92941176 0.92941176 ];
	case 'grey93',
	 colorvector = [ 0.92941176 0.92941176 0.92941176 ];
	case 'gray94',
	 colorvector = [ 0.94117647 0.94117647 0.94117647 ];
	case 'grey94',
	 colorvector = [ 0.94117647 0.94117647 0.94117647 ];
	case 'gray95',
	 colorvector = [ 0.94901961 0.94901961 0.94901961 ];
	case 'grey95',
	 colorvector = [ 0.94901961 0.94901961 0.94901961 ];
	case 'gray96',
	 colorvector = [ 0.96078431 0.96078431 0.96078431 ];
	case 'grey96',
	 colorvector = [ 0.96078431 0.96078431 0.96078431 ];
	case 'gray97',
	 colorvector = [ 0.96862745 0.96862745 0.96862745 ];
	case 'grey97',
	 colorvector = [ 0.96862745 0.96862745 0.96862745 ];
	case 'gray98',
	 colorvector = [ 0.98039216 0.98039216 0.98039216 ];
	case 'grey98',
	 colorvector = [ 0.98039216 0.98039216 0.98039216 ];
	case 'gray99',
	 colorvector = [ 0.98823529 0.98823529 0.98823529 ];
	case 'grey99',
	 colorvector = [ 0.98823529 0.98823529 0.98823529 ];
	case 'gray100',
	 colorvector = [ 1 1 1 ];
	case 'grey100',
	 colorvector = [ 1 1 1 ];
	otherwise
	 colorvector = [ 0.5 0.5 0.5 ];
	 disp('Warning: Color String does not match any known color triple.');
end;

return;
