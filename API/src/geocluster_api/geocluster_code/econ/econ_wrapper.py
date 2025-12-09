from .econ_support import create_teaobject
class Econ():
    def __init__(self, TandP_dict, clgs, drop_down_params, point, econ_params, props_paths, is_heating):
        teaobj = create_teaobject(TandP_dict,
                                        # u_sCO2, u_H2O, c_sCO2, c_H2O,
                                        *clgs,
                                        # case, end_use, fluid, sbt_version,
                                        *drop_down_params,
                                        # mdot, L2, L1, grad, D, Tinj, k,
                                        *point,
                                        # Drilling_cost_per_m, Discount_rate, Lifetime, 
                                        # Direct_use_heat_cost_per_kWth, Power_plant_cost_per_kWe, Pre_Cooling_Delta_T, Turbine_outlet_pressure,
                                        *econ_params,
                                        # properties_H2O_pathname, 
                                        # properties_CO2v2_pathname, 
                                        # additional_properties_CO2v2_pathname,
                                        *props_paths,
                                        is_heating)
        self.tea = teaobj
