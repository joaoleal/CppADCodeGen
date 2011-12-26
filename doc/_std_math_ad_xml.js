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
'std_math_ad.xml'
];
var list_down3 = [
'install.xml',
'introduction.xml',
'ad.xml',
'adfun.xml',
'multi_thread.xml',
'library.xml',
'cppad_ipopt_nlp.xml',
'example.xml',
'preprocessor.xml',
'appendix.xml'
];
var list_down2 = [
'default.xml',
'ad_copy.xml',
'convert.xml',
'advalued.xml',
'boolvalued.xml',
'vecad.xml',
'base_require.xml'
];
var list_down1 = [
'arithmetic.xml',
'std_math_ad.xml',
'mathother.xml',
'condexp.xml',
'discrete.xml',
'user_atomic.xml'
];
var list_down0 = [
'acos.cpp.xml',
'asin.cpp.xml',
'atan.cpp.xml',
'cos.cpp.xml',
'cosh.cpp.xml',
'exp.cpp.xml',
'log.cpp.xml',
'log10.cpp.xml',
'sin.cpp.xml',
'sinh.cpp.xml',
'sqrt.cpp.xml',
'tan.cpp.xml',
'tanh.cpp.xml'
];
var list_current0 = [
'std_math_ad.xml#Syntax',
'std_math_ad.xml#Purpose',
'std_math_ad.xml#x',
'std_math_ad.xml#y',
'std_math_ad.xml#Operation Sequence',
'std_math_ad.xml#fun',
'std_math_ad.xml#Examples',
'std_math_ad.xml#Derivatives',
'std_math_ad.xml#Derivatives.acos',
'std_math_ad.xml#Derivatives.asin',
'std_math_ad.xml#Derivatives.atan',
'std_math_ad.xml#Derivatives.cos',
'std_math_ad.xml#Derivatives.cosh',
'std_math_ad.xml#Derivatives.exp',
'std_math_ad.xml#Derivatives.log',
'std_math_ad.xml#Derivatives.log10',
'std_math_ad.xml#Derivatives.sin',
'std_math_ad.xml#Derivatives.sinh',
'std_math_ad.xml#Derivatives.sqrt',
'std_math_ad.xml#Derivatives.tan',
'std_math_ad.xml#Derivatives.tanh'
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
