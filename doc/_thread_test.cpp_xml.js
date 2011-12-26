var list_across0 = [
'_contents_xml.htm',
'_reference.xml',
'_index.xml',
'_search_xml.htm',
'_external.xml'
];
var list_up0 = [
'cppad.xml',
'multi_thread.xml',
'thread_test.cpp.xml'
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
'parallel_ad.xml',
'team_thread.hpp.xml',
'thread_test.cpp.xml',
'simple_ad.cpp.xml',
'harmonic.cpp.xml',
'multi_newton.cpp.xml'
];
var list_down0 = [
'a11c_openmp.cpp.xml',
'a11c_bthread.cpp.xml',
'a11c_pthread.cpp.xml'
];
var list_current0 = [
'thread_test.cpp.xml#Syntax',
'thread_test.cpp.xml#Running Tests',
'thread_test.cpp.xml#Running Tests.threading',
'thread_test.cpp.xml#Purpose',
'thread_test.cpp.xml#a11c',
'thread_test.cpp.xml#simple_ad',
'thread_test.cpp.xml#harmonic',
'thread_test.cpp.xml#harmonic.test_time',
'thread_test.cpp.xml#harmonic.max_threads',
'thread_test.cpp.xml#harmonic.mega_sum',
'thread_test.cpp.xml#multi_newton',
'thread_test.cpp.xml#multi_newton.test_time',
'thread_test.cpp.xml#multi_newton.max_threads',
'thread_test.cpp.xml#multi_newton.num_zero',
'thread_test.cpp.xml#multi_newton.num_sub',
'thread_test.cpp.xml#multi_newton.num_sum',
'thread_test.cpp.xml#multi_newton.use_ad',
'thread_test.cpp.xml#Threading System Specific Routines',
'thread_test.cpp.xml#Source'
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
