# -*- coding: utf-8 -*-
"""
Created on Sat May  2 13:03:13 2020
@author: d_seet

Integrating a TEV protease-based feedback system for protein aggregation 
burden into the Weiße cell model in a Tellurium reaction-based framework
"""

import tellurium as te
import roadrunner
import antimony 

r = te.loada("""
model Weisse_Aggregate_Burden

/* Tellurium doesn't allow for variables to be entered as reaction
stoichiometry, the values of these parameters (such as n_s, n_tmq) need to be 
entered directly into their reactions */

/* This model uses concentration and time units of # (of molecules) and min */

//Core Weissse Model
    //Reactions
        #Nutrient import & metabolism
        n1: s -> s_i; v_imp;
        n2: s_i -> ; lmda * s_i;
        n3: s_i -> a; v_cat;
        n4: a -> ; lmda * a;
        
        #Transcription & mRNA degradation
        tx1: -> m_r; w_r;
        tx2: m_r -> ; (lmda + d_m) * m_r;
        tx3: -> m_t; w_t;
        tx4: m_t -> ; (lmda + d_m) * m_t;
        tx5: -> m_m; w_m;
        tx6: m_m -> ; (lmda + d_m) * m_m;
        tx7: -> m_q; w_q;
        tx8: m_q -> ; (lmda + d_m) * m_q;
        
        #Ribosome binding
        rb1: m_r + p_r -> c_r; k_b * m_r * p_r;
        rb2: c_r -> m_r + p_r;  k_u * c_r;
        rb3: c_r -> ; lmda * c_r;
        rb4: m_t + p_r -> c_t; k_b * m_t * p_r; 
        rb5: c_t -> m_t + p_r; k_u * c_t;
        rb6: c_t -> ; lmda * c_t;
        rb7: m_m + p_r -> c_m; k_b * m_m * p_r;
        rb8: c_m -> m_m + p_r; k_u * c_m;
        rb9: c_m -> ; lmda * c_m;
        rb10: m_q + p_r -> c_q; k_b * m_q * p_r;
        rb11: c_q -> m_q + p_r; k_u * c_q;
        rb12: c_q -> ; lmda * c_q;
        
        #Translation & protein dilution
        tl1: 7459 a + c_r -> 2 p_r + m_r; v_r;
        tl2: p_r -> ; lmda * p_r;
        tl3: 300 a + c_t -> p_t + m_t + p_r; v_t;
        tl4: p_t -> ; lmda * p_t;
        tl5: 300 a + c_m -> p_m + m_m + p_r; v_m;
        tl6: p_m -> ; lmda * p_m;
        tl7: 300 a + c_q -> p_q + m_q + p_r; v_q;
        tl8: p_q -> ; lmda * p_q;
        
        #Cell population growth
         -> N; lmda * N 
        s -> ; v_imp * (N - 1);  
    
    //Equations
        /* Algebraic relationships and functions are written as Assignement 
        Rules in Antimony. */
    
        #Growth     
        lmda := (gmma / M) * (c_tot);
        gmma := gmma_max * a / (K_gmma + a);
                        
        c_tot := c_r + c_t + c_m + c_q + c_add
        
            /* The term c_add allows for expression of additional proteins to  
            be factored into growth rate. This can be removed if only using
            the base model. */
    
        #Import & metabolism
        v_imp := p_t * (mu_t * s) / (K_t + s);
        v_cat := p_m * (mu_m * s_i) / (K_m + s_i);
        
        #Transcription
        w_r := w_rmax * a / (theta_r + a);
        w_t := w_tmax * a / (theta_nr + a);
        w_m := w_mmax * a / (theta_nr + a);
        w_q := w_qmax * a / (theta_nr + a) / (1 + (p_q / K_q) ^ h_q);
        
        #Translation
        v_r := c_r * gmma / n_r;
        v_t := c_t * gmma / n_tmq;
        v_m := c_m * gmma / n_tmq;
        v_q := c_q * gmma / n_tmq;
    
    //Parameters
        /* Parameter values used are the default values from the 
        Weisse model */
    
        #Growth
        M = 1e8;
        gmma_max = 1260; K_gmma = 7;
        
        #Nutrients
        mu_t = 726; K_t = 1e3; 
        mu_m = 5800; K_m = 1e3; 
        
        #Transcription
        w_rmax = 930; w_tmax = 4.14; w_mmax = 4.14; w_qmax = 949;
        theta_r = 427; theta_nr = 4.38;
        K_q = 152219; h_q = 4; 
        d_m = 0.1;
    
        #Ribosome binding
        k_b = 0.0095; k_u = 1;
        
        #Translation
        n_r = 7459; n_tmq = 300;
        
//Model Extensions
    //Model regulation
        
        /* The 'energy' species s and a occasionally seem to become negative 
        if they approach 0 (due to rounding errors?) and break the model. This
        is a simple fix. */
        
        at (s < 0): s = 0;
        at (a < 0): a = 0;
    
    //Turbidostat
        /* Broadly simulates the behavior of a continuous culture where cell 
        density is kept within a pre-defined range by periodically 
        adding media/removing cells */
        
         -> s; k_in * tur;
        s -> ; d_N * s * tur;
        N -> ; d_N * N * tur;  
        k_in = 5e18; tur = 0;
        d_N := lmda
        
        at (N > 1e10): tur = 1;
        
    //Plotting
        
        /* These aren't part of the model, but are included to facilitate
        plotting graphs */
        
        r_tot := p_r + c_tot;
        r_act := c_tot / r_tot
        p_end := r_tot + p_t + p_m + p_q;
        t_2 := ln(2) / (ln (1 + lmda))

//TEV Feedback Circuit
    /* Includes expression of an inactive TEV protease and a second protein X 
    which forms aggregates. TEV is activated by an input (light) 
    and can then target and cleave X */

    //Reactions    
        #Protein expression
        txa1: -> m_TEV; w_TEV;
        txa2: m_TEV -> ; (lmda + d_m) * m_TEV;
        txa3: -> m_X; w_X * X_induction;
        txa4: m_X -> ; (lmda + d_m) * m_X;
        rba1: m_TEV + p_r -> c_TEV; kb_TEV * m_TEV * p_r;
        rba2: c_TEV -> m_TEV + p_r; ku_TEV * c_TEV;
        rba3: c_TEV -> ; lmda * c_TEV;
        rba4: m_X + p_r -> c_X; kb_X * m_X * p_r;
        rba5: c_X -> m_X + p_r; ku_X * c_X;
        rba6: c_X -> ; lmda * c_X;
        tla1: 242 a + c_TEV -> TEV_in + m_TEV + p_r; v_TEV;
        tla2: TEV_in -> ; lmda * TEV_in;
        tla3: 300 a + c_X -> X_tot + m_X + p_r; v_X;
        tla4: X_tot -> ; lmda * X_tot;
        
        #TEV activation
        uc1: TEV_in -> TEV; k_uc * TEV_in * light_intensity * light_on;
        uc2: TEV -> ; lmda * TEV;
        
        #Protease activity
        pt1: TEV + X_tot -> TX_c; k_on * TEV * X_tot; 
        pt2: TX_c -> TEV + X_tot; k_off * TX_c;
        pt3: TX_c -> ; lmda * TX_c;
        pt4: TX_c -> TEV ; k_cat * TX_c; 
    
    //Equations
        #Protein expression
        w_TEV := w_TEVmax * a / (theta_nr + a);
        w_X := w_Xmax * a / (theta_nr + a);
        v_TEV := c_TEV * gmma / 242; 
        v_X := c_X * gmma / 300;
        
        #Aggregation & Toxicity

        P_Agg := X_tot ^ j / (K_agg ^ j + X_tot ^ j);
        X_agg := X_tot * P_Agg;
        X_sol := X_tot * (1 - P_Agg);
                  
        theta := K_theta ^ h / (K_theta ^ h + X_agg ^ h)
        n3: s_i -> a; v_cat * theta;
        s_i -> ; v_cat * (1 - theta);
            /* As species can't be described by assignement rules, the burden
            a_Eff = theta * a_Weisse needs to be represented by modifying the 
            metabolism reaction. This has the same overall effect. */      
    
    //Parameters  
        #Transcription
        w_TEVmax = 0.8; w_Xmax = 80;
            /* The feedback circuit can be 'removed' from simulations by 
            setting these values to 0 */
        
        #Translation    
        kb_TEV = 0.0095; ku_TEV = 1;
        kb_X = 0.0095; ku_X = 1;
        
        #Protease activity
        k_uc = 1;
        k_on = 2; k_off = 0.1; k_cat = 1.6;
        
        #Aggregate & Toxicity
            /* These parameters probably need to be fit to experimental data.
            They provide a rough indication of aggregation propensity. */
        
        K_agg = 50000; j = 4;
        K_theta = 50000; h = 1; 
        
        #Toggles ('Binary' 0/1 values to switch things on/off)
        X_induction = 0; light_on = 0;
             
    //Events 
        #Start expression of X
            /* This is arbitrarily at t = 300 after reaching the turbidostat
            'steady-state' */
        at (time > 300): X_induction = 1;
        
        #Growth rate monitoring and light activation
            /* This feedback system is activated for 'duration' whenever
            lmda falls below 'threshold' */
        
        at (lmda < threshold): light_on = 1;
        at (timer > duration): light_on = 0;
        at (timer > duration): timer = 0;
         -> timer; light_on;
                
        threshold = 0.018; light_intensity = 1; duration = 5;
            /* Values can be changed to modify the feedback response */
           
    //Plotting
        c_add := c_TEV + c_X
        p_tot := p_end + TEV_in + TEV + X_tot;
        TEV_tot := TEV + TX_c + TEV_in
        TEV_act := TEV + TX_c
        TEV_ratio := (TEV + TX_c) / TEV_tot
        Output := N * X_sol * d_N * tur
        
//Running simulations
    //Cells & nutrients
        s = 1e20; N = 1e7;
    
    //Initializing species
        /* Initial values obtained by running the simulation to steady state 
        without X induction */
        
        m_r = 109; m_t = 16; m_m = 16; m_q = 766;
        c_r = 1105; c_t = 44; c_m = 44; c_q = 2089
        p_r = 1232; p_t = 4423; p_m = 4423; p_q = 2.122e5;
        a = 24;
        
end""")

r.selections = ['time', 'lmda']       #Selects the data to be plotted. 'time' always needs to be present
result = r.simulate(0, 600, 601)        #Simulates model, uses input(start, end, number of points)
te.plotArray(result, xlabel = 'Minutes')
te.show()
