var list_across0 = [
'_contents_xml.htm',
'_reference.xml',
'_index.xml',
'_search_xml.htm',
'_external.xml'
];
var list_up0 = [
'cppad.xml',
'ad.xml',
'advalued.xml',
'arithmetic.xml',
'ad_binary.xml'
];
var list_down3 = [
'default.xml',
'ad_copy.xml',
'convert.xml',
'advalued.xml',
'boolvalued.xml',
'vecad.xml',
'base_require.xml'
];
var list_down2 = [
'arithmetic.xml',
'std_math_ad.xml',
'mathother.xml',
'condexp.xml',
'discrete.xml',
'user_atomic.xml'
];
var list_down1 = [
'unaryplus.xml',
'unaryminus.xml',
'ad_binary.xml',
'compute_assign.xml'
];
var list_down0 = [
'add.cpp.xml',
'sub.cpp.xml',
'mul.cpp.xml',
'div.cpp.xml'
];
var list_current0 = [
'ad_binary.xml#Syntax',
'ad_binary.xml#Purpose',
'ad_binary.xml#Op',
'ad_binary.xml#Base',
'ad_binary.xml#x',
'ad_binary.xml#y',
'ad_binary.xml#z',
'ad_binary.xml#Operation Sequence',
'ad_binary.xml#Example',
'ad_binary.xml#Derivative',
'ad_binary.xml#Derivative.Addition',
'ad_binary.xml#Derivative.Subtraction',
'ad_binary.xml#Derivative.Multiplication',
'ad_binary.xml#Derivative.Division'
];
function choose_across0(item)
{	var index          = item.selectedIndex;
	item.selectedIndex = 0;
	if(index > 0)
		document.location = list_across0[index-1];
}
function choose_up0(item)
{	var index          = item.selectedIndex;
	item.selectedIndex = 0;
	if(index > 0)
		document.location = list_up0[index-1];
}
function choose_down3(item)
{	var index          = item.selectedIndex;
	item.selectedIndex = 0;
	if(index > 0)
		document.location = list_down3[index-1];
}
function choose_down2(item)
{	var index          = item.selectedIndex;
	item.selectedIndex = 0;
	if(index > 0)
		document.location = list_down2[index-1];
}
function choose_down1(item)
{	var index          = item.selectedIndex;
	item.selectedIndex = 0;
	if(index > 0)
		document.location = list_down1[index-1];
}
function choose_down0(item)
{	var index          = item.selectedIndex;
	item.selectedIndex = 0;
	if(index > 0)
		document.location = list_down0[index-1];
}
function choose_current0(item)
{	var index          = item.selectedIndex;
	item.selectedIndex = 0;
	if(index > 0)
		document.location = list_current0[index-1];
}
