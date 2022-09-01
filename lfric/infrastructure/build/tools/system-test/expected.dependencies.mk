# Object dependencies
modules/delta_mod.mod : objects/gamma_mod.o
modules/dir/alpha_mod.mod : objects/dir/alpha_mod.o
modules/dir/zeta_mod.mod : objects/dir/zeta_mod.o
modules/epsilon_mod.mod : objects/epsilon_mod.o
modules/epsilon_one_submod.mod : objects/epsilon_one_submod.o
modules/epsilon_two_submod.mod : objects/epsilon_two_submod.o
modules/eta_mod.mod : objects/eta_mod.o
modules/gamma_mod.mod : objects/gamma_mod.o
modules/mismatched_mod.mod : objects/beta_mod.o
objects/beta_mod.o : modules/delta_mod.mod
objects/dir/alpha_mod.o : modules/eta_mod.mod
objects/epsilon_one_submod.o : modules/epsilon_mod.mod modules/eta_mod.mod
objects/epsilon_two_submod.o : modules/epsilon_mod.mod
objects/one.o : modules/dir/alpha_mod.mod modules/gamma_mod.mod modules/epsilon_mod.mod
objects/two.o : modules/mismatched_mod.mod
