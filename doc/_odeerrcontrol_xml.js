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
'odeerrcontrol.xml'
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
'odeerrcontrol.cpp.xml',
'odeerrmaxabs.cpp.xml'
];
var list_current0 = [
'odeerrcontrol.xml#Syntax',
'odeerrcontrol.xml#Description',
'odeerrcontrol.xml#Include',
'odeerrcontrol.xml#Notation',
'odeerrcontrol.xml#xf',
'odeerrcontrol.xml#Method',
'odeerrcontrol.xml#Method.step',
'odeerrcontrol.xml#Method.Nan',
'odeerrcontrol.xml#Method.order',
'odeerrcontrol.xml#ti',
'odeerrcontrol.xml#tf',
'odeerrcontrol.xml#xi',
'odeerrcontrol.xml#smin',
'odeerrcontrol.xml#smax',
'odeerrcontrol.xml#scur',
'odeerrcontrol.xml#eabs',
'odeerrcontrol.xml#erel',
'odeerrcontrol.xml#ef',
'odeerrcontrol.xml#maxabs',
'odeerrcontrol.xml#nstep',
'odeerrcontrol.xml#Error Criteria Discussion',
'odeerrcontrol.xml#Scalar',
'odeerrcontrol.xml#Vector',
'odeerrcontrol.xml#Example',
'odeerrcontrol.xml#Theory',
'odeerrcontrol.xml#Source Code'
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
