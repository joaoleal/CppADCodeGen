var list_across0 = [
'_contents_xml.htm',
'_reference.xml',
'_index.xml',
'_search_xml.htm',
'_external.xml'
];
var list_up0 = [
'cppad.xml',
'appendix.xml',
'deprecated.xml',
'omp_alloc.xml',
'omp_available.xml'
];
var list_down3 = [
'faq.xml',
'speed.xml',
'theory.xml',
'glossary.xml',
'bib.xml',
'bugs.xml',
'wishlist.xml',
'whats_new.xml',
'deprecated.xml',
'license.xml'
];
var list_down2 = [
'include_deprecated.xml',
'fundeprecated.xml',
'omp_max_thread.xml',
'tracknewdel.xml',
'omp_alloc.xml'
];
var list_down1 = [
'omp_max_num_threads.xml',
'omp_in_parallel.xml',
'omp_get_thread_num.xml',
'omp_get_memory.xml',
'omp_return_memory.xml',
'omp_free_available.xml',
'omp_inuse.xml',
'omp_available.xml',
'omp_create_array.xml',
'omp_delete_array.xml',
'omp_efficient.xml',
'old_max_num_threads.xml',
'omp_alloc.cpp.xml'
];
var list_current0 = [
'omp_available.xml#Deprecated',
'omp_available.xml#Syntax',
'omp_available.xml#Purpose',
'omp_available.xml#thread',
'omp_available.xml#num_bytes',
'omp_available.xml#Example'
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
