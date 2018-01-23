''' --------------------------------------------------------------------------
 *  CppADCodeGen: C++ Algorithmic Differentiation with Source Code Generation:
 *    Copyright (C) 2018 Ciengis
 *
 *  CppADCodeGen is distributed under multiple licenses:
 *
 *   - Eclipse Public License Version 1.0 (EPL1), and
 *   - GNU General Public License Version 3 (GPL3).
 *
 *  EPL1 terms and conditions can be found in the file "epl-v10.txt", while
 *  terms and conditions for the GPL3 can be found in the file "gpl3.txt".
 * ----------------------------------------------------------------------------
 * Author: Joao Leal
'''
import gdb


def is_container(v):
    c = v.type.code
    return c == gdb.TYPE_CODE_STRUCT or c == gdb.TYPE_CODE_UNION


def is_pointer(v):
    return v.type.code == gdb.TYPE_CODE_PTR


class CGPrettyPrinter:
    def __init__(self, val):
        self.val = val

    def to_string(self):
        out = "CG:"

        try:
            node_ptr = self.val['node_']
            if str(node_ptr) != 'NULL':
                node = node_ptr.dereference()
                if str(node['name_']) != '0x0' and str(node['name_']) != 'NULL':
                    out += "\n name: " + str(node['name_'])
                out += "\n variable operation: " + str(node['operation_'])
            else:
                out += "\n <constant>"

            # pointer to value
            pointer = self.val['value_']['_M_t']['_M_t']['_M_head_impl']
            if str(pointer) != 'NULL':
                out += "\n value: " + str(pointer.dereference())

        except Exception as e:
            print('An exception occurred: {}'.format(e))
        except:
            print('An exception occurred')

        return out

    def children(self):
        return [('node', self.val['node_']), ('value', self.val['value_'])]

###############################################################################
#
#  Create and register the pretty printers
#  ---------------------------------------
#
# Register the printers in gdb using the gdb.printing module
#
###############################################################################


def build_pretty_printer():
    pp = gdb.printing.RegexpCollectionPrettyPrinter("CppAD::cg")
    pp.add_printer('CppAD::cg::CG', '^CppAD::cg::CG<.*>$', CGPrettyPrinter)
    return pp


def reload():
    # Remove the pretty printer if it exists
    for printer in gdb.pretty_printers:
        if printer.name == 'CppAD::cg':
            gdb.pretty_printers.remove(printer)
            break

    # Create the new pretty printer
    gdb.printing.register_pretty_printer(gdb.current_objfile(), build_pretty_printer())


reload()
