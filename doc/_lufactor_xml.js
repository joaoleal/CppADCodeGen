var list_across0 = [
'_contents_xml.htm',
'_reference.xml',
'_index.xml',
'_search_xml.htm',
'_external.xml'
];
var list_up0 = [
'cppad.xml',
'library.xml',
'ludetandsolve.xml',
'lufactor.xml'
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
'errorhandler.xml',
'nearequal.xml',
'speed_test.xml',
'speedtest.xml',
'time_test.xml',
'numerictype.xml',
'checknumerictype.xml',
'simplevector.xml',
'checksimplevector.xml',
'nan.xml',
'pow_int.xml',
'poly.xml',
'ludetandsolve.xml',
'rombergone.xml',
'rombergmul.xml',
'runge45.xml',
'rosen34.xml',
'odeerrcontrol.xml',
'odegear.xml',
'odegearcontrol.xml',
'benderquad.xml',
'opt_val_hes.xml',
'luratio.xml',
'cppad_vector.xml',
'thread_alloc.xml',
'memory_leak.xml'
];
var list_down1 = [
'lusolve.xml',
'lufactor.xml',
'luinvert.xml'
];
var list_down0 = [
'lufactor.cpp.xml',
'lu_factor.hpp.xml'
];
var list_current0 = [
'lufactor.xml#Syntax',
'lufactor.xml#Description',
'lufactor.xml#Include',
'lufactor.xml#Matrix Storage',
'lufactor.xml#sign',
'lufactor.xml#ip',
'lufactor.xml#jp',
'lufactor.xml#LU',
'lufactor.xml#LU.A',
'lufactor.xml#LU.P',
'lufactor.xml#LU.L',
'lufactor.xml#LU.U',
'lufactor.xml#LU.Factor',
'lufactor.xml#LU.Determinant',
'lufactor.xml#SizeVector',
'lufactor.xml#FloatVector',
'lufactor.xml#Float',
'lufactor.xml#AbsGeq',
'lufactor.xml#Example',
'lufactor.xml#Source'
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
