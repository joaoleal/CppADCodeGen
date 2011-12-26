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
'user_atomic.xml',
'mat_mul.cpp.xml',
'mat_mul.hpp.xml'
];
var list_down3 = [
'arithmetic.xml',
'std_math_ad.xml',
'mathother.xml',
'condexp.xml',
'discrete.xml',
'user_atomic.xml'
];
var list_down2 = [
'user_tan.cpp.xml',
'mat_mul.cpp.xml'
];
var list_down1 = [
'mat_mul.hpp.xml'
];
var list_current0 = [
'mat_mul.hpp.xml#Syntax',
'mat_mul.hpp.xml#Example',
'mat_mul.hpp.xml#Begin Source',
'mat_mul.hpp.xml#Extra Call Information',
'mat_mul.hpp.xml#Matrix Indexing',
'mat_mul.hpp.xml#One Matrix Multiply',
'mat_mul.hpp.xml#Reverse Partials One Order',
'mat_mul.hpp.xml#Set Union',
'mat_mul.hpp.xml#CppAD User Atomic Callback Functions',
'mat_mul.hpp.xml#Declare mat_mul Function'
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
