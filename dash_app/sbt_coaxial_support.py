def check_coaxial_diameters(Diameter1, Diameter2, PipeParam3, A_flow_min):

    # PipeParam3 = center pipe thickness | thickness = float(PipeParam3)

    # Guardrail: coaxial geometry requires inner (center pipe) diameter < outer (wellbore/annulus) diameter.
    # Also ensure that outerradiuscenterpipe (radiuscenterpipe + thicknesscenterpipe) < radius
    # If UI / stored values violate this, SBT2 can end up with negative/invalid annulus hydraulic diameter,
    # which shows up as negative Reynolds numbers and "laminar flow" termination.
    # TODO AB: An SME needs to evaluate this.
    if Diameter1 is not None and Diameter2 is not None:
        try:
            d1 = float(Diameter1) # Diameter1 = wellbore
            d2 = float(Diameter2)  # Diameter2 = center_pipe
            thickness_val = float(PipeParam3)
            
            # Calculate radii (after set_tube_geometry divides by 2)
            radius_wellbore = d1 / 2
            radiuscenterpipe = d2 / 2
            outerradiuscenterpipe = radiuscenterpipe + thickness_val
            
            if outerradiuscenterpipe >= radius_wellbore:
                """ 
                Need to reduce either Diameter2 or thicknesscenterpipe. 
                Strategy: Ensure minimum flow area for annulus (not just clearance). Minimum annulus flow area: 1e-4 m² = 100 cm² .
                For annulus: A = π*(r_wellbore² - r_outer²) . Solving: r_outer_max = sqrt(r_wellbore² - A_min/π)
                """
                r_outer_max_sq = radius_wellbore**2 - A_flow_min / np.pi
                if r_outer_max_sq <= 0:
                    # Wellbore is too small to accommodate minimum flow area
                    error_msg = (f"Error: Wellbore diameter ({d1:.6f} m) is too small to accommodate "
                               f"minimum annulus flow area ({A_flow_min:.6e} m²). "
                               f"Minimum wellbore diameter required: {2*np.sqrt(A_flow_min/np.pi):.6f} m. "
                               f"Simulation terminated.")
                    # print(f"[ERROR] {error_msg}", flush=True)
                    raise ValueError(error_msg)
                max_outerradius = np.sqrt(r_outer_max_sq)
                # Also ensure at least 5mm clearance as a safety margin
                min_annulus_clearance = 0.005  # 5mm minimum clearance
                max_outerradius = min(max_outerradius, radius_wellbore - min_annulus_clearance)
                
                if thickness_val > max_outerradius:
                    # Thickness is too large, clamp it
                    new_thickness = max(0.005, max_outerradius - radiuscenterpipe)
                    if new_thickness < 0.005:
                        # Even with minimum thickness, we need to reduce center pipe radius
                        new_radiuscenterpipe = max_outerradius - 0.005
                        new_d2 = new_radiuscenterpipe * 2
                        Diameter2 = new_d2
                    else:
                        PipeParam3 = new_thickness
                else:
                    # Thickness is OK, but center pipe radius is too large
                    new_radiuscenterpipe = max_outerradius - thickness_val
                    new_d2 = new_radiuscenterpipe * 2
                    Diameter2 = new_d2
            
            # Also check if Diameter2 >= Diameter1 (before accounting for thickness)
            elif d2 >= d1:
                # Use a more reasonable clamp: 40% of wellbore diameter, but at least 0.15m
                new_d2 = max(0.15, min(0.4 * d1, d1 * 0.4))
                Diameter2 = new_d2
            
            # After any clamping, validate that flow areas are reasonable
            # Recalculate to get final values
            final_radius_wellbore = float(Diameter1) / 2
            final_radiuscenterpipe = float(Diameter2) / 2
            final_thickness = float(PipeParam3)
            final_outerradius = final_radiuscenterpipe + final_thickness
            
            # Calculate flow areas
            A_flow_annulus = np.pi * (final_radius_wellbore**2 - final_outerradius**2)
            A_flow_centerpipe = np.pi * final_radiuscenterpipe**2
            
            # Ensure minimum flow area to prevent extremely high velocities
            # CO2 requires larger flow areas due to lower density (~200-800 kg/m³ vs ~1000 kg/m³ for H2O)
            # For CO2, we need ~5x larger flow area to keep velocities reasonable
            if self.fluid == "sCO2":
                A_flow_min = 5e-3  # Minimum flow area for CO2 [m²] = 5000 cm² (5x larger than H2O)
                A_flow_min_reason = "CO2 has much lower density than H2O (~200-800 kg/m³ vs ~1000 kg/m³), requiring larger flow areas to prevent excessive velocities"
            else:
                A_flow_min = 1e-4  # Minimum flow area for H2O [m²] = 100 cm²
                A_flow_min_reason = "standard minimum for H2O"
            
            if A_flow_annulus < A_flow_min or A_flow_centerpipe < A_flow_min:
                error_msg = (f"Error: After geometry clamping, flow areas are too small for numerical stability ({self.fluid}). "
                           f"A_flow_annulus={A_flow_annulus:.6e} m², A_flow_centerpipe={A_flow_centerpipe:.6e} m². "
                           f"Minimum required: {A_flow_min:.6e} m² ({A_flow_min_reason}). "
                           f"Diameter1={Diameter1:.6f} m, Diameter2={Diameter2:.6f} m, "
                           f"thickness={final_thickness:.6f} m. "
                           f"This geometry combination is not suitable for {self.fluid} simulation. "
                           f"Consider: (1) increasing wellbore diameter, (2) reducing center pipe diameter, or (3) reducing mass flow rate. "
                           f"Simulation terminated.")
                print(f"[ERROR] {error_msg}", flush=True)
                raise ValueError(error_msg)
            
            # Pre-check: Estimate velocity for CO2 to warn about potential issues
            # CO2 has much lower density than H2O, so velocities will be higher
            if self.fluid == "sCO2" and mdot is not None:
                # Estimate minimum CO2 density (gas phase can be ~1-10 kg/m³, supercritical ~200-800 kg/m³)
                # Use conservative estimate of 200 kg/m³ for supercritical CO2, but also check gas phase (1 kg/m³)
                min_co2_density_supercritical = 200.0  # kg/m³ (conservative estimate for supercritical CO2)
                min_co2_density_gas = 1.0  # kg/m³ (worst case: gas phase CO2)
                
                # Check both annulus and center pipe velocities
                estimated_vel_annulus_supercritical = mdot / (A_flow_annulus * min_co2_density_supercritical)
                estimated_vel_annulus_gas = mdot / (A_flow_annulus * min_co2_density_gas)
                estimated_vel_centerpipe_supercritical = mdot / (A_flow_centerpipe * min_co2_density_supercritical)
                estimated_vel_centerpipe_gas = mdot / (A_flow_centerpipe * min_co2_density_gas)
                
                max_estimated_velocity = max(estimated_vel_annulus_gas, estimated_vel_centerpipe_gas)
                velocity_error_threshold = 600.0  # m/s - error if estimated velocity exceeds this
                velocity_warning_threshold = 400.0  # m/s - warn if estimated velocity exceeds this
                
                # If gas phase velocity exceeds limit, this geometry won't work
                if max_estimated_velocity > velocity_error_threshold:
                    error_msg = (f"Error: Estimated CO2 velocity too high for this geometry. "
                               f"Estimated max velocity (gas phase): {max_estimated_velocity:.2f} m/s > {velocity_error_threshold:.2f} m/s limit. "
                               f"A_flow_annulus={A_flow_annulus:.6e} m², A_flow_centerpipe={A_flow_centerpipe:.6e} m², mdot={mdot} kg/s. "
                               f"CO2 density can be very low (~1-10 kg/m³ in gas phase), causing extremely high velocities. "
                               f"Consider: (1) increasing wellbore diameter, (2) reducing center pipe diameter, or (3) reducing mass flow rate. "
                               f"Simulation terminated.")
                    print(f"[ERROR] {error_msg}", flush=True)
                    raise ValueError(error_msg)
                
                if max_estimated_velocity > velocity_warning_threshold:
                    warning_msg = (f"[WARNING] Pre-check: Estimated CO2 velocity may be high (max={max_estimated_velocity:.1f} m/s) "
                                 f"for geometry with A_flow_annulus={A_flow_annulus:.6e} m², A_flow_centerpipe={A_flow_centerpipe:.6e} m², mdot={mdot} kg/s. "
                                 f"CO2 density (~200-800 kg/m³ supercritical, ~1-10 kg/m³ gas phase) is much lower than H2O (~1000 kg/m³), causing higher velocities. "
                                 f"Annulus (gas phase): {estimated_vel_annulus_gas:.1f} m/s, Center pipe (gas phase): {estimated_vel_centerpipe_gas:.1f} m/s. "
                                 f"Consider increasing flow area or reducing mass flow rate to avoid numerical instability.")
                    print(warning_msg, flush=True)
        except (TypeError, ValueError) as e:
            # If parsing fails, let SBT handle validation downstream.
            print(f"[WARNING] Could not validate coaxial geometry: {e}", flush=True)
            pass
    

# code from sbt_utils that could integrate around line 314:

        # # Validate geometry to prevent numerical instability
        # if A_flow_annulus <= 0:
        #     raise ValueError(f"Invalid annulus flow area: {A_flow_annulus:.6f} m². "
        #                    f"radius={radius:.6f} m, outerradiuscenterpipe={outerradiuscenterpipe:.6f} m. "
        #                    f"This indicates invalid coaxial geometry (center pipe outer radius >= wellbore radius).")
        # if Dh_annulus <= 0:
        #     raise ValueError(f"Invalid annulus hydraulic diameter: {Dh_annulus:.6f} m. "
        #                    f"radius={radius:.6f} m, outerradiuscenterpipe={outerradiuscenterpipe:.6f} m. "
        #                    f"This indicates invalid coaxial geometry (center pipe outer radius >= wellbore radius).")
        # if A_flow_centerpipe <= 0:
        #     raise ValueError(f"Invalid center pipe flow area: {A_flow_centerpipe:.6f} m². "
        #                    f"radiuscenterpipe={radiuscenterpipe:.6f} m. "
        #                    f"This indicates invalid center pipe geometry.")
        # # Warn if annulus is very thin (may cause numerical instability)
        # min_annulus_area_ratio = 0.01  # 1% of wellbore area
        # min_annulus_area = min_annulus_area_ratio * math.pi * radius**2
        # if A_flow_annulus < min_annulus_area:
        #     print(f"[WARNING] Annulus flow area is very small ({A_flow_annulus:.6f} m² < {min_annulus_area:.6f} m²). "
        #           f"This may cause numerical instability in the solver. "
        #           f"radius={radius:.6f} m, outerradiuscenterpipe={outerradiuscenterpipe:.6f} m", flush=True)
        