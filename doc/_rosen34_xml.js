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
'rosen34.xml'
];
var list_down2 = [
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
var list_down1 = [
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
var list_down0 = [
'rosen34.cpp.xml'
];
var list_current0 = [
'rosen34.xml#Syntax',
'rosen34.xml#Description',
'rosen34.xml#Include',
'rosen34.xml#xf',
'rosen34.xml#Fun',
'rosen34.xml#Fun.t',
'rosen34.xml#Fun.x',
'rosen34.xml#Fun.f',
'rosen34.xml#Fun.f_t',
'rosen34.xml#Fun.f_x',
'rosen34.xml#Fun.Nan',
'rosen34.xml#Fun.Warning',
'rosen34.xml#Fun.Optimization',
'rosen34.xml#M',
'rosen34.xml#ti',
'rosen34.xml#tf',
'rosen34.xml#xi',
'rosen34.xml#e',
'rosen34.xml#Scalar',
'rosen34.xml#Vector',
'rosen34.xml#Parallel Mode',
'rosen34.xml#Example',
'rosen34.xml#Source Code'
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
