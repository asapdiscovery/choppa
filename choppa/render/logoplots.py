import pandas as pd
import logomaker
import matplotlib.pyplot as plt
import matplotlib
import base64
import io 
import numpy as np
import math

LOGOPLOT_WITHOUT_CONF_COLORSCHEME = { # see https://github.com/jbkinney/logomaker/blob/master/logomaker/src/colors.py
        'A': '#f76ab4',
        'C': '#ff7f00',
        'D': '#e41a1c',
        'E': '#e41a1c',
        'F': '#84380b',
        'G': '#f76ab4',
        'H': '#3c58e5',
        'I': '#12ab0d',
        'K': '#3c58e5',
        'L': '#12ab0d',
        'M': '#12ab0d',
        'N': '#972aa8',
        'P': '#12ab0d',
        'Q': '#972aa8',
        'R': '#3c58e5',
        'S': '#ff7f00',
        'T': '#ff7f00',
        'V': '#12ab0d',
        'W': '#84380b',
        'Y': '#84380b',
        'X': '#000000', # add 'X' so that LogoMaker doesn't log to stdout
    }
WHITE_EMPTY_SQUARE = "iVBORw0KGgoAAAANSUhEUgAAAJYAAACfCAIAAACUbLd9AAAACXBIWXMAAAsTAAALEwEAmpwYAAABhElEQVR4nO3RwQkAIBDAMHX/nc8hfEghmaDQPTOLsvM7gFcW5lmYZ2GehXkW5lmYZ2GehXkW5lmYZ2GehXkW5lmYZ2GehXkW5lmYZ2GehXkW5lmYZ2GehXkW5lmYZ2GehXkW5lmYZ2GehXkW5lmYZ2GehXkW5lmYZ2GehXkW5lmYZ2GehXkW5lmYZ2GehXkW5lmYZ2GehXkW5lmYZ2GehXkW5lmYZ2GehXkW5lmYZ2GehXkW5lmYZ2GehXkW5lmYZ2GehXkW5lmYZ2GehXkW5lmYZ2GehXkW5lmYZ2GehXkW5lmYZ2GehXkW5lmYZ2GehXkW5lmYZ2GehXkW5lmYZ2GehXkW5lmYZ2GehXkW5lmYZ2GehXkW5lmYZ2GehXkW5lmYZ2GehXkW5lmYZ2GehXkW5lmYZ2GehXkW5lmYZ2GehXkW5lmYZ2GehXkW5lmYZ2GehXkW5lmYZ2GehXkW5lmYZ2GehXkW5lmYZ2GehXkW5lmYZ2GehXkW5lmYZ2GehXkW5l1WGwQ7i50I0AAAAABJRU5ErkJggg=="

def render_singleres_logoplot(res):
    """
    Renders a simple, single-letter 'logoplot', normally for showing the wildtype residue. The plot
    is square and is typically rendered top center downstream.
    """
    _, ax = plt.subplots(figsize=(4, 4))

    # create Logo object
    logomaker.Logo(
        pd.DataFrame({res:1}, index=[0]),
        font_name="Sans Serif",
        color_scheme=LOGOPLOT_WITHOUT_CONF_COLORSCHEME,
        flip_below=False,
        show_spines=True,
        ax=ax
    )

    plt.xticks([])
    plt.yticks([])
    # plt.savefig("debug_logoplot.png", dpi=70, bbox_inches="tight") # uncomment for testing
    # plt object directly to base64 string instead of tmpfile
    lp_bytes = io.BytesIO()
    plt.savefig(
        lp_bytes,  
        format='png', 
        dpi=70, # DPI 70 seems to be ~smallest we can get away with
        bbox_inches="tight",
        )
    lp_bytes.seek(0)
    lp_base64 = base64.b64encode(lp_bytes.read())
    plt.close()

    return lp_base64 

class LogoPlot():
    """
    Given a dict with mutants for a given residue, generate logoplots.
    """
    def __init__(self, residue_dict, fitness_threshold):
        self.residue_dict = residue_dict
        self.fitness_threshold = fitness_threshold

    def divide_fitness_types(self):
        """
        Determines which mutants are fit/unfit given the `fitness_threshold` and returns 
        all data required for logoplot generation in a simple dict form.
        """
        # first collect some info on the wildtype residue
        wildtype = self.residue_dict['wildtype']['aa']
        wildtype_fitness = [self.residue_dict['wildtype']['fitness']] # make this a list to keep consistent typing

        if not math.isnan(self.residue_dict['wildtype']['confidence']):
            # if the wildtype has a confidence tag we can assume that all mutant data will also have confidence
            self.confidence = True
            wildtype_fitness.append(self.residue_dict['wildtype']['confidence'])            
        else:
            self.confidence = False

        # get fit and unfit mutants according to `fitness_threshold`
        unfit_mutants = {}
        fit_mutants = {}
        for mutant in self.residue_dict['mutants']:
            if self.confidence: # also write confidence into the dict
                if mutant['fitness'] < self.fitness_threshold: 
                    unfit_mutants[mutant['aa']] = [mutant['fitness'], mutant['confidence']]
                else:
                    fit_mutants[mutant['aa']] = [mutant['fitness'], mutant['confidence']]
            elif not self.confidence: # no confidences, so only write fitness to the dict
                if mutant['fitness'] < self.fitness_threshold: 
                    unfit_mutants[mutant['aa']] = [mutant['fitness']] # make this a list to keep consistent typing
                else:
                    fit_mutants[mutant['aa']] = [mutant['fitness']]

        return {wildtype: wildtype_fitness}, unfit_mutants, fit_mutants  

    def render_logoplot(self, mutants, global_min_confidence=False, global_max_confidence=False, lhs=True, wildtype=False):
        """
        Creates a logoplot as a base64 string. Also annotes with confidence values if present.

        TODO: nicer rounded ticks agnostic to array limits
        """  
        if len(mutants) == 0:
            # this can happen when there are no mutants in this category. Return an empty white-sqare base64 instead.
            return WHITE_EMPTY_SQUARE
        plt.switch_backend('Agg') # prevents plt from opening a figure on OS
        if wildtype: # we want this to be a bit smaller and square because it'll always have 1 residue.
            _, ax = plt.subplots(figsize=(4, 4))
        else:
            _, ax = plt.subplots(figsize=(3, 10))

        # if there are confidences, we well color the logoplot AA letters by confidence and
        # show a color bar if this is the left-hand-side logoplot. 
        if self.confidence:
            if not global_min_confidence or not global_max_confidence:
                raise ValueError("If confidence is provided then a global confidence limit needs to be passed to render_logoplot()")
            matplotlib.rcParams.update({'font.size': 12}) # instead of doing for ticks/title separately

            # define a 'mappable' which allows us to generate the colormap before the rest of the plot
            mappable = plt.cm.ScalarMappable(
                cmap=matplotlib.colors.LinearSegmentedColormap.from_list(
                    'custom blue', # bit convoluted but this way we force the colormap to be continuous
                    ['#ff6600','#0066ff'], # between two colors
                    N=256
                    ), 
                norm=matplotlib.colors.Normalize(vmin=global_min_confidence, vmax=global_max_confidence),
                )
            if not lhs:
                # plot the colorbar
                plt.colorbar(
                    mappable, 
                    ticks=np.linspace(global_min_confidence, global_max_confidence, int(len(mutants)/4)), 
                    ax=ax,
                ).ax.set_title("         Confidence", y=1.02) # indent to make title appear nicely centered. `ha` doesn't get us there

            # build a dict that has {residue : RGBA color, ..} that we can use to color the logoplot
            conf_color_per_AA = { k:mappable.to_rgba(v[1]) for k,v in mutants.items() }
        else:
            # just use regular coloring if there is no confidence set.
            conf_color_per_AA = LOGOPLOT_WITHOUT_CONF_COLORSCHEME

        # create Logo object
        logomaker.Logo(
            pd.DataFrame(mutants)[:1],
            font_name="Sans Serif",
            color_scheme=conf_color_per_AA,
            flip_below=False,
            show_spines=True,
            ax=ax
        )

        plt.xticks([])
        plt.yticks([])
        # plt.savefig("debug_logoplot.png", dpi=70, bbox_inches="tight") # uncomment for testing
        # plt object directly to base64 string instead of tmpfile
        lp_bytes = io.BytesIO()
        plt.savefig(
            lp_bytes,  
            format='png', 
            dpi=70, # DPI 70 seems to be ~smallest we can get away with
            bbox_inches="tight",
            )
        lp_bytes.seek(0)
        lp_base64 = base64.b64encode(lp_bytes.read())
        plt.close()

        return lp_base64

    def build_logoplot(self, global_min_confidence=False, global_max_confidence=False):
        # determine the wildtype, unfit and fit mutants for this input
        wildtype, unfit_mutants, fit_mutants = self.divide_fitness_types()

        # generate the logoplot base64 for wildtype (LHS, top), fit (LHS, bottom) and unfit (RHS; with colorbar) 
        wildtype_base64 = self.render_logoplot(wildtype, global_min_confidence=global_min_confidence, 
                             global_max_confidence=global_max_confidence, wildtype=True)
        fit_base64 = self.render_logoplot(fit_mutants, global_min_confidence=global_min_confidence, 
                             global_max_confidence=global_max_confidence)
        unfit_base64 = self.render_logoplot(unfit_mutants, global_min_confidence=global_min_confidence, 
                             global_max_confidence=global_max_confidence, lhs=False)
        
        return wildtype_base64, fit_base64, unfit_base64
        
if __name__ == "__main__":
    # test a fitness dict with conf values
    residue_dict = {'fitness_aligned_index': 164, 'fitness_csv_index': 160, 'wildtype': {'aa': 'L', 'fitness': 1.0, 'confidence': 4221}, 'mutants': [{'aa': 'V', 'fitness': -1.98, 'confidence': 2455}, {'aa': 'I', 'fitness': -3.3, 'confidence': 434}, {'aa': 'E', 'fitness': -4.52, 'confidence': 3706}, {'aa': 'Q', 'fitness': -3.78, 'confidence': 3079}, {'aa': 'D', 'fitness': 0.56, 'confidence': 3615}, {'aa': 'N', 'fitness': -1.05, 'confidence': 3911}, {'aa': 'H', 'fitness': -0.59, 'confidence': 4891}, {'aa': 'W', 'fitness': -1.88, 'confidence': 2627}, {'aa': 'F', 'fitness': -0.56, 'confidence': 2663}, {'aa': 'Y', 'fitness': 0.66, 'confidence': 4534}, {'aa': 'R', 'fitness': -0.73, 'confidence': 11}, {'aa': 'K', 'fitness': 0.89, 'confidence': 3763}, {'aa': 'S', 'fitness': -1.77, 'confidence': 2352}, {'aa': 'T', 'fitness': -1.16, 'confidence': 3843}, {'aa': 'M', 'fitness': -2.11, 'confidence': 4018}, {'aa': 'A', 'fitness': -3.33, 'confidence': 4132}, {'aa': 'G', 'fitness': -0.44, 'confidence': 3251}, {'aa': 'P', 'fitness': -0.0, 'confidence': 3817}, {'aa': 'C', 'fitness': -1.5, 'confidence': 3445}, {'aa': 'X', 'fitness': -0.37, 'confidence': 3281}]}
    LogoPlot(residue_dict, fitness_threshold=0.7).build_logoplot(global_min_confidence=10, global_max_confidence=5000)

    # test a fitness dict without conf values
    # nan=np.nan
    # residue_dict = {'fitness_aligned_index': 165, 'fitness_csv_index': 161, 'wildtype': {'aa': 'V', 'fitness': 1.0, 'confidence': nan}, 'mutants': [{'aa': 'I', 'fitness': 0.09, 'confidence': nan}, {'aa': 'L', 'fitness': -2.13, 'confidence': nan}, {'aa': 'E', 'fitness': -2.3, 'confidence': nan}, {'aa': 'Q', 'fitness': -1.96, 'confidence': nan}, {'aa': 'D', 'fitness': -3.96, 'confidence': nan}, {'aa': 'N', 'fitness': -2.51, 'confidence': nan}, {'aa': 'H', 'fitness': -3.09, 'confidence': nan}, {'aa': 'W', 'fitness': -3.94, 'confidence': nan}, {'aa': 'F', 'fitness': -2.76, 'confidence': nan}, {'aa': 'Y', 'fitness': -1.46, 'confidence': nan}, {'aa': 'R', 'fitness': -4.67, 'confidence': nan}, {'aa': 'K', 'fitness': -3.52, 'confidence': nan}, {'aa': 'S', 'fitness': -0.39, 'confidence': nan}, {'aa': 'T', 'fitness': -4.71, 'confidence': nan}, {'aa': 'M', 'fitness': -4.5, 'confidence': nan}, {'aa': 'A', 'fitness': -0.8, 'confidence': nan}, {'aa': 'G', 'fitness': -3.9, 'confidence': nan}, {'aa': 'P', 'fitness': 0.4, 'confidence': nan}, {'aa': 'C', 'fitness': -4.86, 'confidence': nan}, {'aa': 'X', 'fitness': -1.61, 'confidence': nan}]}
    # print(LogoPlot(residue_dict, fitness_threshold=0.7).build_logoplot(global_min_confidence=10, global_max_confidence=5000))

