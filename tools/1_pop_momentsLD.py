import moments.LD

with open("params_1_pop.csv") as fin:
    params = fin.readlines()
  
p_names = params[0].strip().split(",")

with open("moments_LD_1_pop_out.csv", "w+") as fout:
    stats = moments.LD.Util.moment_names(1)
    l = ",".join(stats[0] + stats[1]) + "\n"
    fout.write(l)
    for p_list in params[1:]:
        p_list = p_list.split(",")
        N = float(p_list[p_names.index('N')])
        u = float(p_list[p_names.index('mu')])
        r = float(p_list[p_names.index('r')])
        y = moments.LD.Demographics1D.snm(rho = 4 * N * r, theta = 4 * N * u)
        l = ",".join([str(_) for _ in y.LD()[0]]) + "\n"
        fout.write(l)
