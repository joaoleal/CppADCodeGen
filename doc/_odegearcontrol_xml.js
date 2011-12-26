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
'odegearcontrol.xml'
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
'odegearcontrol.cpp.xml'
];
var list_current0 = [
'odegearcontrol.xml#Syntax',
'odegearcontrol.xml#Purpose',
'odegearcontrol.xml#Include',
'odegearcontrol.xml#Notation',
'odegearcontrol.xml#xf',
'odegearcontrol.xml#Fun',
'odegearcontrol.xml#Fun.t',
'odegearcontrol.xml#Fun.x',
'odegearcontrol.xml#Fun.f',
'odegearcontrol.xml#Fun.f_x',
'odegearcontrol.xml#Fun.Warning',
'odegearcontrol.xml#M',
'odegearcontrol.xml#ti',
'odegearcontrol.xml#tf',
'odegearcontrol.xml#xi',
'odegearcontrol.xml#smin',
'odegearcontrol.xml#smax',
'odegearcontrol.xml#sini',
'odegearcontrol.xml#eabs',
'odegearcontrol.xml#erel',
'odegearcontrol.xml#ef',
'odegearcontrol.xml#maxabs',
'odegearcontrol.xml#nstep',
'odegearcontrol.xml#Error Criteria Discussion',
'odegearcontrol.xml#Scalar',
'odegearcontrol.xml#Vector',
'odegearcontrol.xml#Example',
'odegearcontrol.xml#Theory',
'odegearcontrol.xml#Source Code'
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
