import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def plot_blowoff_analysis(filename='test_data.csv'):
    # Load data
    df = pd.read_csv(filename)
    
    # Set style
    sns.set(style="whitegrid")
    
    # --- PLOT 1: Design Viability (Da vs Pressure Drop) ---
    plt.figure(figsize=(10, 6))
    sns.scatterplot(data=df, x='dP_percent', y='Da', hue='Stabilization', style='Stabilization', alpha=0.7)
    
    # Add Constraints
    plt.axvline(x=2.0, color='r', linestyle='--', label='Max Pressure Drop (2%)')
    plt.axhline(y=1.0, color='k', linestyle='--', label='Blowoff Limit (Approx Da=1)')
    
    plt.title('Design Viability Map: Stability vs. Efficiency')
    plt.xlabel('Pressure Drop (%)')
    plt.ylabel('Damköhler Number (Da)')
    plt.legend()
    plt.grid(True)
    plt.savefig('plot_viability.png')
    plt.show()

    # --- PLOT 2: Operating Window (Da vs Phi) ---
    plt.figure(figsize=(10, 6))
    # Filter for a specific viable U0 (e.g., 100 m/s) to make the plot clearer
    subset = df[df['U0'] == df['U0'].median()] 
    
    sns.lineplot(data=subset, x='PHI', y='Da', hue='Stabilization')
    
    plt.axhline(y=1.0, color='r', linestyle='--', label='Unstable Region')
    plt.title(f'Stability Operating Window (at U0 ~ {subset.U0.iloc[0]:.0f} m/s)')
    plt.xlabel('Equivalence Ratio ($\phi$)')
    plt.ylabel('Stability Margin (Da)')
    plt.legend()
    plt.savefig('plot_operating_window.png')
    plt.show()

    # --- PLOT 3: Jet Flame Speed Check ---
    # Filter for Jet Flames only
    jet_df = df[df['Stabilization'] == 'Jet Flame'].copy()
    
    if not jet_df.empty:
        plt.figure(figsize=(10, 6))
        
        # Calculate absolute Turbulent Flame Speed (ST) 
        # ST_SL_* columns are ratios, so multiply by S_L
        jet_df['ST_Zimont_Abs'] = jet_df['ST_SL_Zimont'] * jet_df['S_L']
        
        sns.scatterplot(data=jet_df, x='U0', y='ST_Zimont_Abs', hue='PHI')
        
        # Add a 1:1 line (Theoretical holding limit where Flow = Flame Speed)
        # Note: Actual blowout happens when Flow >> Flame Speed, but this is a reference
        max_val = max(jet_df['U0'].max(), jet_df['ST_Zimont_Abs'].max())
        plt.plot([0, max_val], [0, max_val], 'k--', label='Velocity Balance (U = ST)')
        
        plt.title('Jet Flame: Aerodynamic Stretch vs. Flame Speed')
        plt.xlabel('Flow Velocity U0 (m/s)')
        plt.ylabel('Turbulent Flame Speed ST (m/s)')
        plt.legend(title='Phi')
        plt.savefig('plot_jet_stability.png')
        plt.show()

if __name__ == "__main__":
    try:
        plot_blowoff_analysis("/Users/raja/Library/CloudStorage/OneDrive-SharedLibraries-McGillUniversity/CombustionDynamicsF2025 - General/04. Code - git repo/src/test_data.csv")
    except FileNotFoundError:
        print("Error: 'test_data.csv' not found. Run your main simulation first.")