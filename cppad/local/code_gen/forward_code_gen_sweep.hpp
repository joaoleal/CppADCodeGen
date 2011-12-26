/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2011 Ciengis

CppAD is distributed under multiple licenses. This distribution is under
the terms of the 
                    Common Public License Version 1.0.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */
#ifndef CPPAD_CODE_GEN_FORWARD_CODE_GEN_SWEEP_INCLUDED
#define CPPAD_CODE_GEN_FORWARD_CODE_GEN_SWEEP_INCLUDED

CPPAD_BEGIN_NAMESPACE

#define CPPAD_FORWARD_SWEEP_TRACE 0

template <class Base>
size_t forward_code_gen_sweep(
std::ostream& s_out,
CodeGenNameProvider<Base>& names,
bool print,
size_t d,
size_t n,
size_t numvar,
player<Base> *Rec) {
    // op code for current instruction
    OpCode op;

    // index for current instruction
    size_t i_op;

    // next variables 
    size_t i_var;

#if CPPAD_USE_FORWARD0SWEEP
    CPPAD_ASSERT_UNKNOWN(d > 0);
#else
    addr_t* non_const_arg;
#endif
    const addr_t* arg = 0;

    // temporary indices
    size_t i, ell;

    // initialize the comparison operator (ComOp) counter
    size_t compareCount = 0;

    // if this is an order zero calculation, initialize vector indices
    pod_vector<size_t> VectorInd; // address for each element
    pod_vector<bool> VectorVar; // is element a variable
    i = Rec->num_rec_vecad_ind();
    if (i > 0) {
        VectorInd.extend(i);
        VectorVar.extend(i);
        while (i--) {
            VectorInd[i] = Rec->GetVecInd(i);
            VectorVar[i] = false;
        }
    }

    // work space used by UserOp.
    const size_t user_k = d; // order of this forward mode calculation
    const size_t user_k1 = d + 1; // number of orders for this calculation
    vector<Base> user_tx; // argument vector Taylor coefficients 
    vector<Base> user_ty; // result vector Taylor coefficients 
    vector<size_t> user_iy; // variable indices for results vector
    size_t user_index = 0; // indentifier for this user_atomic operation
    size_t user_id = 0; // user identifier for this call to operator
    size_t user_i = 0; // index in result vector
    size_t user_j = 0; // index in argument vector
    size_t user_m = 0; // size of result vector
    size_t user_n = 0; // size of arugment vector
    // next expected operator in a UserOp sequence

    enum {
        user_start, user_arg, user_ret, user_end
    } user_state = user_start;

    // check numvar argument
    CPPAD_ASSERT_UNKNOWN(Rec->num_rec_var() == numvar);

    // length of the parameter vector (used by CppAD assert macros)
    const size_t num_par = Rec->num_rec_par();

    // pointer to the beginning of the parameter vector
    const Base* parameter = 0;
    if (num_par > 0)
        parameter = Rec->GetPar();

#if ! CPPAD_USE_FORWARD0SWEEP
    // length of the text vector (used by CppAD assert macros)
    const size_t num_text = Rec->num_rec_text();

    // pointer to the beginning of the text vector
    const char* text = 0;
    if (num_text > 0)
        text = Rec->GetTxt(0);
#endif

    // skip the BeginOp at the beginning of the recording
    Rec->start_forward(op, arg, i_op, i_var);
    CPPAD_ASSERT_UNKNOWN(op == BeginOp);
#if CPPAD_FORWARD_SWEEP_TRACE
    std::cout << std::endl;
#endif
    while (op != EndOp) {
        // this op
        Rec->next_forward(op, arg, i_op, i_var);
        CPPAD_ASSERT_UNKNOWN((i_op > n) | (op == InvOp));
        CPPAD_ASSERT_UNKNOWN((i_op <= n) | (op != InvOp));

        // action depends on the operator
        switch (op) {
            case AbsOp:
                forward_code_gen_abs_op(s_out, names, d, i_var, arg[0]);
                break;
                // -------------------------------------------------

            case AddvvOp:
                forward_code_gen_addvv_op(s_out, names, d, i_var, arg, parameter);
                break;
                // -------------------------------------------------

            case AddpvOp:
                CPPAD_ASSERT_UNKNOWN(size_t(arg[0]) < num_par);
                forward_code_gen_addpv_op(s_out, names, d, i_var, arg, parameter);
                break;
                // -------------------------------------------------

            case AcosOp:
                // sqrt(1 - x * x), acos(x)
                CPPAD_ASSERT_UNKNOWN(i_var < numvar);
                forward_code_gen_acos_op(s_out, names, d, i_var, arg[0]);
                break;
                // -------------------------------------------------

            case AsinOp:
                // sqrt(1 - x * x), asin(x)
                CPPAD_ASSERT_UNKNOWN(i_var < numvar);
                forward_code_gen_asin_op(s_out, names, d, i_var, arg[0]);
                break;
                // -------------------------------------------------

            case AtanOp:
                // 1 + x * x, atan(x)
                CPPAD_ASSERT_UNKNOWN(i_var < numvar);
                forward_code_gen_atan_op(s_out, names, d, i_var, arg[0]);
                break;
                // -------------------------------------------------

                //			case CSumOp:
                //			// CSumOp has a variable number of arguments and
                //			// next_forward thinks it one has one argument.
                //			// we must inform next_forward of this special case.
                //			Rec->forward_csum(op, arg, i_op, i_var);
                //			forward_code_gen_csum_op(
                //				d, i_var, arg, num_par, parameter, J, Taylor
                //			);
                //			break;
                // -------------------------------------------------

            case CExpOp:
                forward_code_gen_cond_op(s_out, names, d, i_var, arg, num_par, parameter);
                break;
                // ---------------------------------------------------

            case ComOp:
#if ! USE_FORWARD0SWEEP 
                if (d == 0)
                    forward_code_gen_comp_op_0(s_out, names, arg, num_par, parameter);
#endif
                break;
                // ---------------------------------------------------

            case CosOp:
                // sin(x), cos(x)
                CPPAD_ASSERT_UNKNOWN(i_var < numvar);
                forward_code_gen_cos_op(s_out, names, d, i_var, arg[0]);
                break;
                // ---------------------------------------------------

            case CoshOp:
                // sinh(x), cosh(x)
                CPPAD_ASSERT_UNKNOWN(i_var < numvar);
                forward_code_gen_cosh_op(s_out, names, d, i_var, arg[0]);
                break;
                // -------------------------------------------------

            case DisOp:
                forward_code_gen_dis_op<Base > (s_out, names, d, i_var, arg);
                break;
                // -------------------------------------------------

            case DivvvOp:
                forward_code_gen_divvv_op(s_out, names, d, i_var, arg);
                break;
                // -------------------------------------------------

            case DivpvOp:
                CPPAD_ASSERT_UNKNOWN(size_t(arg[0]) < num_par);
                forward_code_gen_divpv_op(s_out, names, d, i_var, arg, parameter);
                break;
                // -------------------------------------------------

            case DivvpOp:
                CPPAD_ASSERT_UNKNOWN(size_t(arg[1]) < num_par);
                forward_code_gen_divvp_op(s_out, names, d, i_var, arg, parameter);
                break;
                // -------------------------------------------------

            case EndOp:
                CPPAD_ASSERT_UNKNOWN(NumArg(op) == 0);
                CPPAD_ASSERT_UNKNOWN(NumRes(op) == 0);
                //CPPAD_ASSERT_NARG_NRES(op, 0, 0);
                break;
                // -------------------------------------------------

            case ExpOp:
                forward_code_gen_exp_op(s_out, names, d, i_var, arg[0]);
                break;
                // -------------------------------------------------
                //
            case InvOp:
                CPPAD_ASSERT_UNKNOWN(NumArg(op) == 0);
                break;
                // -------------------------------------------------
                //
                //			case LdpOp:
                //# if ! CPPAD_USE_FORWARD0SWEEP
                //			if( d == 0 )
                //			{	non_const_arg = Rec->forward_non_const_arg();
                //				forward_load_p_op_0(
                //					i_var, 
                //					non_const_arg, 
                //					num_par, 
                //					parameter, 
                //					J, 
                //					Taylor,
                //					Rec->num_rec_vecad_ind(),
                //					VectorVar.data(),
                //					VectorInd.data()
                //				);
                //			}
                //			else
                //# endif
                //			{	forward_load_op( op, d, i_var, arg, J, Taylor);
                //			}
                //			break;
                //			// -------------------------------------------------
                //
                //			case LdvOp:
                //# if ! CPPAD_USE_FORWARD0SWEEP
                //			if( d == 0 )
                //			{	non_const_arg = Rec->forward_non_const_arg();
                //				forward_load_v_op_0(
                //					i_var, 
                //					non_const_arg, 
                //					num_par, 
                //					parameter, 
                //					J, 
                //					Taylor,
                //					Rec->num_rec_vecad_ind(),
                //					VectorVar.data(),
                //					VectorInd.data()
                //				);
                //			}
                //			else
                //# endif
                //			{	forward_code_gen_load_op( op, d, i_var, arg, J, Taylor);
                //			}
                //			break;
                // -------------------------------------------------

            case LogOp:
                forward_code_gen_log_op(s_out, names, d, i_var, arg[0]);
                break;
                // -------------------------------------------------

            case MulvvOp:
                forward_code_gen_mulvv_op(s_out, names, d, i_var, arg);
                break;
                // -------------------------------------------------

            case MulpvOp:
                CPPAD_ASSERT_UNKNOWN(size_t(arg[0]) < num_par);
                forward_code_gen_mulpv_op(s_out, names, d, i_var, arg, parameter);
                break;
                // -------------------------------------------------

            case ParOp:
                forward_code_gen_par_op(s_out, names, d, i_var, arg, num_par, parameter);
                break;
                // -------------------------------------------------

            case PowvpOp:
                CPPAD_ASSERT_UNKNOWN(size_t(arg[1]) < num_par);
                forward_code_gen_powvp_op(s_out, names, d, i_var, arg, parameter);
                break;
                // -------------------------------------------------

            case PowpvOp:
                CPPAD_ASSERT_UNKNOWN(size_t(arg[0]) < num_par);
                forward_code_gen_powpv_op(s_out, names, d, i_var, arg, parameter);
                break;
                // -------------------------------------------------

            case PowvvOp:
                forward_code_gen_powvv_op(s_out, names, d, i_var, arg);
                break;
                // -------------------------------------------------

            case PriOp:
#if ! CPPAD_USE_FORWARD0SWEEP
                if (print) forward_code_gen_pri_0(s_out, names,
                        i_var, arg, num_text, text, num_par, parameter);
#endif
                break;
                // -------------------------------------------------

            case SinOp:
                // cos(x), sin(x)
                CPPAD_ASSERT_UNKNOWN(i_var < numvar);
                forward_code_gen_sin_op(s_out, names, d, i_var, arg[0]);
                break;
                // -------------------------------------------------

            case SinhOp:
                // cosh(x), sinh(x)
                CPPAD_ASSERT_UNKNOWN(i_var < numvar);
                forward_code_gen_sinh_op(s_out, names, d, i_var, arg[0]);
                break;
                // -------------------------------------------------

            case SqrtOp:
                forward_code_gen_sqrt_op(s_out, names, d, i_var, arg[0]);
                break;
                // -------------------------------------------------
                //
                //			case StppOp:
                //# if ! CPPAD_USE_FORWARD0SWEEP
                //			if( d == 0 )
                //			{	forward_code_gen_store_pp_op_0(
                //					i_var, 
                //					arg, 
                //					num_par, 
                //					J, 
                //					Taylor,
                //					Rec->num_rec_vecad_ind(),
                //					VectorVar.data(),
                //					VectorInd.data()
                //				);
                //			}
                //# endif
                //			break;
                //			// -------------------------------------------------
                //
                //			case StpvOp:
                //# if ! CPPAD_USE_FORWARD0SWEEP
                //			if( d == 0 )
                //			{	forward_code_gen_store_pv_op_0(
                //					i_var, 
                //					arg, 
                //					num_par, 
                //					J, 
                //					Taylor,
                //					Rec->num_rec_vecad_ind(),
                //					VectorVar.data(),
                //					VectorInd.data()
                //				);
                //			}
                //# endif
                //			break;
                //			// -------------------------------------------------
                //
                //			case StvpOp:
                //# if ! CPPAD_USE_FORWARD0SWEEP
                //			if( d == 0 )
                //			{	forward_code_gen_store_vp_op_0(
                //					i_var, 
                //					arg, 
                //					num_par, 
                //					J, 
                //					Taylor,
                //					Rec->num_rec_vecad_ind(),
                //					VectorVar.data(),
                //					VectorInd.data()
                //				);
                //			}
                //# endif
                //			break;
                //			// -------------------------------------------------
                //
                //			case StvvOp:
                //# if ! CPPAD_USE_FORWARD0SWEEP
                //			if( d == 0 )
                //			{	forward_code_gen_store_vv_op_0(
                //					i_var, 
                //					arg, 
                //					num_par, 
                //					J, 
                //					Taylor,
                //					Rec->num_rec_vecad_ind(),
                //					VectorVar.data(),
                //					VectorInd.data()
                //				);
                //			}
                //# endif
                //			break;
                // -------------------------------------------------

            case SubvvOp:
                forward_code_gen_subvv_op(s_out, names, d, i_var, arg, parameter);
                break;
                // -------------------------------------------------

            case SubpvOp:
                CPPAD_ASSERT_UNKNOWN(size_t(arg[0]) < num_par);
                forward_code_gen_subpv_op(s_out, names, d, i_var, arg, parameter);
                break;
                // -------------------------------------------------

            case SubvpOp:
                CPPAD_ASSERT_UNKNOWN(size_t(arg[1]) < num_par);
                forward_code_gen_subvp_op(s_out, names, d, i_var, arg, parameter);
                break;
                // -------------------------------------------------

            case TanOp:
                // tan(x)^2, tan(x)
                CPPAD_ASSERT_UNKNOWN(i_var < numvar);
                forward_code_gen_tan_op(s_out, names, d, i_var, arg[0]);
                break;
                // -------------------------------------------------

            case TanhOp:
                // tanh(x)^2, tanh(x)
                CPPAD_ASSERT_UNKNOWN(i_var < numvar);
                forward_code_gen_tanh_op(s_out, names, d, i_var, arg[0]);
                break;
                // -------------------------------------------------

                //			case UserOp:
                //			// start or end an atomic operation sequence
                //			CPPAD_ASSERT_UNKNOWN( NumRes( UserOp ) == 0 );
                //			CPPAD_ASSERT_UNKNOWN( NumArg( UserOp ) == 4 );
                //			if( user_state == user_start )
                //			{	user_index = arg[0];
                //				user_id    = arg[1];
                //				user_n     = arg[2];
                //				user_m     = arg[3];
                //				if(user_tx.size() < user_n * user_k1)
                //					user_tx.resize(user_n * user_k1);
                //				if(user_ty.size() < user_m * user_k1)
                //					user_ty.resize(user_m * user_k1);
                //				if(user_iy.size() < user_m)
                //					user_iy.resize(user_m);
                //				user_j     = 0;
                //				user_i     = 0;
                //				user_state = user_arg;
                //			}
                //			else
                //			{	CPPAD_ASSERT_UNKNOWN( user_state == user_end );
                //				CPPAD_ASSERT_UNKNOWN( user_index == size_t(arg[0]) );
                //				CPPAD_ASSERT_UNKNOWN( user_id    == size_t(arg[1]) );
                //				CPPAD_ASSERT_UNKNOWN( user_n     == size_t(arg[2]) );
                //				CPPAD_ASSERT_UNKNOWN( user_m     == size_t(arg[3]) );
                //				user_state = user_start;
                //
                //				// call users function for this operation
                //				user_atomic<Base>::forward(user_index, user_id,
                //					user_k, user_n, user_m, user_tx, user_ty
                //				);
                //				for(i = 0; i < user_m; i++) if( user_iy[i] > 0 )
                //					Taylor[ user_iy[i] * J + user_k ] = 
                //						user_ty[ i * user_k1 + user_k ];
                //			}
                //			break;
                //
                //			case UsrapOp:
                //			// parameter argument in an atomic operation sequence
                //			CPPAD_ASSERT_UNKNOWN( user_state == user_arg );
                //			CPPAD_ASSERT_UNKNOWN( user_j < user_n );
                //			CPPAD_ASSERT_UNKNOWN( size_t(arg[0]) < num_par );
                //			user_tx[user_j * user_k1 + 0] = parameter[ arg[0]];
                //			for(ell = 1; ell < user_k1; ell++)
                //				user_tx[user_j * user_k1 + ell] = Base(0);
                //			++user_j;
                //			if( user_j == user_n )
                //				user_state = user_ret;
                //			break;
                //
                //			case UsravOp:
                //			// variable argument in an atomic operation sequence
                //			CPPAD_ASSERT_UNKNOWN( user_state == user_arg );
                //			CPPAD_ASSERT_UNKNOWN( user_j < user_n );
                //			CPPAD_ASSERT_UNKNOWN( size_t(arg[0]) <= i_var );
                //			for(ell = 0; ell < user_k1; ell++)
                //				user_tx[user_j * user_k1 + ell] = Taylor[ arg[0] * J + ell];
                //			++user_j;
                //			if( user_j == user_n )
                //				user_state = user_ret;
                //			break;
                //
                //			case UsrrpOp:
                //			// parameter result in an atomic operation sequence
                //			CPPAD_ASSERT_UNKNOWN( user_state == user_ret );
                //			CPPAD_ASSERT_UNKNOWN( user_i < user_m );
                //			user_iy[user_i] = 0;
                //			user_ty[user_i * user_k1 + 0] = parameter[ arg[0]];
                //			for(ell = 1; ell < user_k; ell++)
                //				user_ty[user_i * user_k1 + ell] = Base(0);
                //			user_i++;
                //			if( user_i == user_m )
                //				user_state = user_end;
                //			break;
                //
                //			case UsrrvOp:
                //			// variable result in an atomic operation sequence
                //			CPPAD_ASSERT_UNKNOWN( user_state == user_ret );
                //			CPPAD_ASSERT_UNKNOWN( user_i < user_m );
                //			user_iy[user_i] = i_var;
                //			for(ell = 0; ell < user_k; ell++)
                //				user_ty[user_i * user_k1 + ell] = Taylor[ i_var * J + ell];
                //			user_i++;
                //			if( user_i == user_m )
                //				user_state = user_end;
                //			break;
                // -------------------------------------------------

            default:
                CPPAD_ASSERT_UNKNOWN(0);
        }
#if CPPAD_FORWARD_SWEEP_TRACE
        size_t i_tmp = i_var;
        Base* Z_tmp = Taylor + J * i_var;
        printOp(
                std::cout,
                Rec,
                i_tmp,
                op,
                arg,
                d + 1,
                Z_tmp,
                0,
                (Base *) CPPAD_NULL
                );
    }
    std::cout << std::endl;
#else
    }
#endif
    CPPAD_ASSERT_UNKNOWN(user_state == user_start);
    CPPAD_ASSERT_UNKNOWN(i_var + 1 == Rec->num_rec_var());

    return compareCount;
}

#undef CPPAD_FORWARD_SWEEP_TRACE

CPPAD_END_NAMESPACE
#endif
