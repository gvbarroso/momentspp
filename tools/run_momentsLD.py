import moments.LD

with open("params_2_pops.csv") as fin:
    params = fin.readlines()

p_names = params[0].strip().split(",")

with open("moments_LD_out.csv", "w+") as fout:
    stats = moments.LD.Util.moment_names(2)
    l = ",".join(stats[0] + stats[1]) + "\n"
    fout.write(l)
    for p_list in params[1:]:
        p_list = p_list.split(",")
        u = eval(p_list[p_names.index('mu')])
        r = eval(p_list[p_names.index('r')])
        m01 = eval(p_list[p_names.index('m_01')])
        m10 = eval(p_list[p_names.index('m_10')])
        N0 = eval(p_list[p_names.index('N_0')])
        N1 = eval(p_list[p_names.index('N_1')])
        rho = 4 * N0 * r
        theta = 4 * N0 * u
        nus = [N0 / N0, N1 / N0]
        m = [[0, 2 * N0 * m01], [2 * N0 * m10, 0]]
        y = moments.LD.LDstats(
            moments.LD.Numerics.steady_state_two_pop(nus, m, rho=rho, theta=theta),
            num_pops=2,
        )
        num_epochs = eval(p_list[p_names.index('num_epochs')])
        if num_epochs == 2:
            T = eval(p_list[p_names.index('num_gen_epoch_2')]) / 2 / N1
            N0b = eval(p_list[p_names.index('N_0_epoch_2')])
            N1b = eval(p_list[p_names.index('N_1_epoch_2')])
            nus = [N0b / N0, N1b / N0]
            y.integrate(nus, T, theta=theta, rho=rho, m=m)
        else:
            assert num_epochs == 1
        l = ",".join([str(_) for _ in y[0]] + [str(_) for _ in y[1]]) + "\n"
        fout.write(l)
