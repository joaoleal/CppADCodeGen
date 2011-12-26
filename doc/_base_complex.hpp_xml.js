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
'base_require.xml',
'base_complex.hpp.xml'
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
'base_cond_exp.xml',
'base_identical.xml',
'base_ordered.xml',
'base_std_math.xml',
'base_float.hpp.xml',
'base_double.hpp.xml',
'base_complex.hpp.xml',
'base_adolc.hpp.xml'
];
var list_down0 = [
'complexpoly.cpp.xml',
'not_complex_ad.cpp.xml'
];
var list_current0 = [
'base_complex.hpp.xml#Example',
'base_complex.hpp.xml#See Also',
'base_complex.hpp.xml#Include Order',
'base_complex.hpp.xml#CondExpOp',
'base_complex.hpp.xml#CondExpRel',
'base_complex.hpp.xml#EqualOpSeq',
'base_complex.hpp.xml#Identical',
'base_complex.hpp.xml#Ordered',
'base_complex.hpp.xml#Integer',
'base_complex.hpp.xml#epsilon',
'base_complex.hpp.xml#isnan',
'base_complex.hpp.xml#Unary Standard Math',
'base_complex.hpp.xml#pow'
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
