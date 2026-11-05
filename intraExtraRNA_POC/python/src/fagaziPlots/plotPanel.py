"""
Andrew E. Davidson
aedavids@ucsc.edu
Date: 2025-09-26
Copyright (c) 2020-2023, Regents of the University of California All rights reserved.
https://polyformproject.org/licenses/noncommercial/1.0.0

Description: 

    functions to make creating manuscript final artwork complient figures
    
    def checkForFonts(
                    fontList : list[str] = ["Arial", "Helvetica", "DejaVu Sans", "Aptos"]
                ) -> None: 
    
    def combineFigures(
            figures : list[plt.figure],
            outputStem : str ="combinedFigure",
            widthMm=180, 
            heightMm=185, 
            nRows=3,
            nCols=2,
            addPanelLabels=True,
            bottomMargin=None,
        ) -> list[Path, Path]:
        
    @dataclass
    class Panel:
        label: str
        leftMm: float
        topMm: float
        widthMm: float
        heightMm: float
        figure: Figure | None = None
        
    def combineFiguresLayout(
        panels : dict[str, Panel],
        outputStem="combinedFigure",
        widthMm=180,
        heightMm=134,
        addPanelLabels=True,
        labelFontSize=7,
        legendFigure=None,
        legendHeightMm=10,
        pngDpi=600,
    ) -> list[Path, Path] :
    
    def displayPDFInNotebook(
        pdfPath : Path,
        dpi : int = 150, 
        width=900
        ) -> None :
        
    def getSubplotWidthAndHeight(layout :dict[str, Panel],
                             panelId : str ) -> tuple[float, float]:
                             
    def removeLegends(
        listOfFigs : list[plt.figure],
        panelWidthInMM : int,
        legendHeightMm : int,
        nCols : int,
    ) -> tuple[list[plt.figure], plt.figure]: 
    
    example of use src/jupyterNotebooks/binary.vs.multclassifiers.ipynb
    
    natureManuscriptStyle = {
        "font.family": "Arial",
        "font.size": 7,
        "axes.titlesize": 7,
        "axes.labelsize": 7,
        "xtick.labelsize": 7,
        "ytick.labelsize": 7,
        "legend.fontsize": 7,
        "legend.title_fontsize": 7,
        "figure.titlesize": 7,
        "pdf.fonttype": 42,
        "savefig.bbox": None,
    }
    
    def setFigureFont(fig, fontFamily="Arial"):
            ''' 
            make sure all figures use the correct final artwork compliant font
            
            example of use
            
               manuscriptStyle = {
                        "font.family": "Arial",
                        "font.size": 7,
                        "axes.titlesize": 7,
                        "axes.labelsize": 7,
                        "xtick.labelsize": 7,
                        "ytick.labelsize": 7,
                        "legend.fontsize": 7,
                        "legend.title_fontsize": 7,
                        "figure.titlesize": 7,
                        "pdf.fonttype": 42,
                        "savefig.bbox": None,
                    }
                    
            
            with mpl.rc_context(manuscriptStyle):
                fig, ax = yourPlottingFunction(...)
                setFigureFont(fig)
                # Apply your existing layout adjustment here.
            '''     

    see deprecated make_panel()
        arguments a list containing ConfusionMatrixDisplay, Figures, or axes object
        return a single figure.
    
    
    
    
generate by chatgpt
"""

import logging
logger = logging.getLogger(__name__)

from copy import deepcopy
from dataclasses import dataclass
import fitz

import io
from io import BytesIO
from itertools import combinations
# from IPython.display import Image, display

import math
from math import isfinite
from typing import List, Optional, Union, Literal, Sequence

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from matplotlib import font_manager
from matplotlib.figure import Figure
from matplotlib.text import Text

from pathlib import Path
import pymupdf
from sklearn.metrics import ConfusionMatrixDisplay
from types import MappingProxyType # create read only mappint

# https://www.nature.com/documents/NRJs-guide-to-preparing-final-artwork.pdf
natureManuscriptStyle = MappingProxyType({
                            # 10/15/26 not sure what happend Arial is not avalible on mustard
                            # DejaVu Sans is matplotlib default. should be complient
                            # missing font generates a ton of warnigns
                            "font.family": "DejaVu Sans", #"Arial",
                            "font.size": 7,
                            "axes.titlesize": 7,
                            "axes.labelsize": 7,
                            "xtick.labelsize": 7,
                            "ytick.labelsize": 7,
                            "legend.fontsize": 7,
                            "legend.title_fontsize": 7,
                            "figure.titlesize": 7,
                            "pdf.fonttype": 42,
                            "savefig.bbox": None,
                        })
   

@dataclass
class Panel:
    label: str
    leftMm: float
    topMm: float
    widthMm: float
    heightMm: float
    figure: Figure | None = None
    
################################################################################
def checkForFonts(
        fontList : list[str] = ["Arial", "Helvetica", "DejaVu Sans", "Aptos"]
    ) -> None: 
    '''
    check to see that matplotlib fonts are installed on local machine
    '''
    for fontName in fontList:
        try:
            fontPath = font_manager.findfont(
                font_manager.FontProperties(family=fontName),
                fallback_to_default=False,
            )
            print(f"{fontName}: {fontPath}")
        except ValueError:
            print(f"{fontName}: not available")

################################################################################
def combineFigures(
    figures : list[plt.figure],
    outputStem : str ="combinedFigure",
    widthMm=180, 
    heightMm=185, 
    nRows=3,
    nCols=2,
    addPanelLabels=True,
    bottomMargin=None,
    legendFigure=None,
    legendHeightMm=None,
) -> list[Path, Path]:
    """
    Assemble Matplotlib figures in row order into PDF and PNG files.
    
    TODO: create a version that works well when grid cell size vary
    
    Changing the default grid can change the displayed font sizes, 
    because each source figure is scaled to fit its position. 
    To preserve the original font sizes, create each input figure at:
    
    figsize = (
        (widthMm / nCols) / 25.4,
        (heightMm / nRows) / 25.4,
    )
    
    arguments 
        outputStem :
            example f'{runConfig.imgPanelOutDir}/foo/myPanel                 
            will save a version with suffix '.png' and '.pdf'
        
        widthMm=180, nature biomedical engineering 2 text column width
        
        heightMm=185, Chosen figure height; accommodates our three-row layout
        The artwork guide gives approximately 185 mm as the maximum two-column 
        figure height with a caption under 300 words.
        
        https://www.nature.com/documents/NRJs-guide-to-preparing-final-artwork.pdf
        
        bottomMargin :
                use to adjust the bottom margin x-axis label is cliped
                example bottomMargin=0.20
                
        legendFigure=None,
            an optional fig to place on the bottom of the panel

        legendHeightMm=18,
            create some extra space along the bottom for a panel level legend

    returns (pdfPath, pngPath)
    
    Optional legend example see commented out function at bottom of file
    def exampleOfCombinedFiguresWithOptionalLegend()
            
    """
    figures = list(figures)

    if nRows < 1 or nCols < 1:
        raise ValueError("nRows and nCols must be positive integers.")
    
    if widthMm <= 0 or heightMm <= 0:
        raise ValueError("Figure width and height must be positive.")
    
    if not figures:
        raise ValueError("Provide at least one figure.")
    
    if len(figures) > nRows * nCols:
        raise ValueError("The grid has fewer positions than figures.")
    
    if bottomMargin is not None and not 0 <= bottomMargin < 1:
        raise ValueError("bottomMargin must be between 0 and 1.")
    
    if legendFigure is not None and legendHeightMm <= 0:
        raise ValueError("legendHeightMm must be positive.")

    outputStem = Path(outputStem).expanduser()
    outputStem.parent.mkdir(parents=True, exist_ok=True)

    pdfPath = Path(f"{outputStem}.pdf")
    pngPath = Path(f"{outputStem}.png")

    mmToPoints = 72 / 25.4
    pageWidth = widthMm * mmToPoints
    
    # heightMm describes only the data grid.
    gridHeight = heightMm * mmToPoints
    
    legendHeight = (
        legendHeightMm * mmToPoints
        if legendFigure is not None
        else 0
    )
    
    pageHeight = gridHeight + legendHeight

    panelWidth = pageWidth / nCols
    panelHeight = gridHeight / nRows
    
    # pageHeight = heightMm * mmToPoints
    # panelWidth = pageWidth / nCols
    # panelHeight = pageHeight / nRows

    def panelLetter(index):
        # a–z, then aa, ab, ...
        label = ""
        number = index + 1
        while number:
            number, remainder = divmod(number - 1, 26)
            label = chr(ord("a") + remainder) + label
        return label

    with pymupdf.open() as combined:
        page = combined.new_page(
            width=pageWidth,
            height=pageHeight,
        )

        for index, fig in enumerate(figures):
            row, column = divmod(index, nCols)

            rectangle = pymupdf.Rect(
                column * panelWidth,
                row * panelHeight,
                (column + 1) * panelWidth,
                (row + 1) * panelHeight,
            )

            originalBottom = fig.subplotpars.bottom

            with BytesIO() as buffer:
                try:
                    if bottomMargin is not None:
                        fig.subplots_adjust(bottom=bottomMargin)

                    with mpl.rc_context({
                        "savefig.bbox": None,
                        "pdf.fonttype": 42,
                    }):
                        fig.savefig(
                            buffer,
                            format="pdf",
                            bbox_inches=None,
                            facecolor="white",
                        )
                finally:
                    if bottomMargin is not None:
                        fig.subplots_adjust(bottom=originalBottom)

                with pymupdf.open(
                    stream=buffer.getvalue(),
                    filetype="pdf",
                ) as source:
                    page.show_pdf_page(
                        rectangle,
                        source,
                        0,
                        keep_proportion=True,
                    )

            if addPanelLabels:
                page.insert_text(
                    (
                        rectangle.x0 + 2 * mmToPoints,
                        rectangle.y0 + 4 * mmToPoints,
                    ),
                    panelLetter(index),
                    fontname="hebo",  # Helvetica Bold
                    fontsize=7,
                    color=(0, 0, 0),
                )

        

        # Place the shared legend below the data grid.
        if legendFigure is not None:
            legendRectangle = pymupdf.Rect(
                0,
                gridHeight,
                pageWidth,
                pageHeight,
            )

            with BytesIO() as buffer:
                with mpl.rc_context({
                    "savefig.bbox": None,
                    "pdf.fonttype": 42,
                }):
                    legendFigure.savefig(
                        buffer,
                        format="pdf",
                        bbox_inches=None,
                        facecolor="white",
                    )

                with pymupdf.open(
                    stream=buffer.getvalue(),
                    filetype="pdf",
                ) as source:
                    page.show_pdf_page(
                        legendRectangle,
                        source,
                        0,
                        keep_proportion=True,
                    )
        
        combined.save(str(pdfPath))
        page.get_pixmap(dpi=600, alpha=False).save(str(pngPath))

    return pdfPath, pngPath


################################################################################
def combineFiguresLayout(
        panels : dict[str, Panel],
        outputStem="combinedFigure",
        widthMm=180,
        heightMm=134,
        addPanelLabels=True,
        labelFontSize=7,
        legendFigure=None,
        legendHeightMm=10,
        pngDpi=600,
        maxWidthMm=180,
        maxHeightMm=185,
    ) ->  tuple[Path, Path]:
    """
    Assemble figures using explicit rectangles in millimetres.
    
    argument 
        nat bio med eng final artwork max for 2 col text
        maxWidthMm default = 180
        maxHeightMm default = 185



    Each panel is a dictionary containing:
        figure: Matplotlib Figure
        rectMm: (left, top, width, height)
        label: optional panel label, e.g. "a"
        
        example
        panels = {
            "a": Panel("a", leftMm=0, topMm=0, widthMm=60, heightMm=65, figure=f1),
            "b": Panel("d", leftMm=14.5, topMm=71.5, widthMm=80, heightMm=60, figure=f2),
        }

    Coordinates start at the top-left of the data area.
    heightMm excludes the optional legend strip.

    Source figures should match their rectangle dimensions.
    Mismatches produce a warning; placement preserves proportions.

    Returns:
        (pdfPath, pngPath)
    """
    validatePanelLayout(
                    panels=panels,
                    widthMm=widthMm,
                    heightMm=heightMm,
                    legendFigure=legendFigure,
                    legendHeightMm=legendHeightMm,
                    maxWidthMm=maxWidthMm,
                    maxHeightMm=maxHeightMm,
                    )
    # Convert the keyed layout into the internal list used below.
    panelList = []

    for panel in panels.values():
        if panel.figure is None:
            raise ValueError(f"Panel {panel.label!r} has no figure.")

        panelList.append({
            "label": panel.label,
            "figure": panel.figure,
            "rectMm": (
                panel.leftMm,
                panel.topMm,
                panel.widthMm,
                panel.heightMm,
            ),
        })
        
    panels = panelList

    if not panels:
        raise ValueError("Provide at least one panel.")
    if widthMm <= 0 or heightMm <= 0:
        raise ValueError("Figure dimensions must be positive.")
    if legendFigure is not None and legendHeightMm <= 0:
        raise ValueError("legendHeightMm must be positive.")

    # Validate placements before writing any files.
    for index, panel in enumerate(panels):
        left, top, width, height = panel["rectMm"]

        if width <= 0 or height <= 0:
            raise ValueError(f"Panel {index}: dimensions must be positive.")

        if (
            left < 0
            or top < 0
            or left + width > widthMm + 1e-6
            or top + height > heightMm + 1e-6
        ):
            raise ValueError(f"Panel {index} extends outside the data area.")

    outputStem = Path(outputStem).expanduser()
    outputStem.parent.mkdir(parents=True, exist_ok=True)

    pdfPath = Path(f"{outputStem}.pdf")
    pngPath = Path(f"{outputStem}.png")

    mmToPoints = 72 / 25.4
    totalHeightMm = heightMm + (
        legendHeightMm if legendFigure is not None else 0
    )

    def placeFigure(page, figure, rectMm, name):
        left, top, width, height = rectMm
        sourceWidth, sourceHeight = figure.get_size_inches() * 25.4

        if (
            abs(sourceWidth - width) > 0.1
            or abs(sourceHeight - height) > 0.1
        ):
            logger.warnings.warn(
                f"{name}: source is "
                f"{sourceWidth:.1f} × {sourceHeight:.1f} mm; "
                f"destination is {width:.1f} × {height:.1f} mm. "
                "Placement will scale the content, including fonts.",
                stacklevel=2,
            )

        rectangle = pymupdf.Rect(
            left * mmToPoints,
            top * mmToPoints,
            (left + width) * mmToPoints,
            (top + height) * mmToPoints,
        )

        with BytesIO() as buffer:
            with mpl.rc_context({
                "savefig.bbox": None,
                "pdf.fonttype": 42,
            }):
                figure.savefig(
                    buffer,
                    format="pdf",
                    bbox_inches=None,
                    facecolor="white",
                )

            with pymupdf.open(
                stream=buffer.getvalue(),
                filetype="pdf",
            ) as source:
                page.show_pdf_page(
                    rectangle,
                    source,
                    0,
                    keep_proportion=True,
                )

    with pymupdf.open() as combined:
        page = combined.new_page(
            width=widthMm * mmToPoints,
            height=totalHeightMm * mmToPoints,
        )

        for index, panel in enumerate(panels):
            placeFigure(
                page,
                panel["figure"],
                panel["rectMm"],
                f"Panel {index + 1}",
            )

            label = panel.get("label")
            if addPanelLabels and label:
                left, top, _, _ = panel["rectMm"]

                page.insert_text(
                    (
                        (left + 2) * mmToPoints,
                        (top + 4) * mmToPoints,
                    ),
                    label,
                    fontname="hebo",  # Helvetica Bold
                    fontsize=labelFontSize,
                    color=(0, 0, 0),
                )

        if legendFigure is not None:
            placeFigure(
                page,
                legendFigure,
                (0, heightMm, widthMm, legendHeightMm),
                "Shared legend",
            )

        combined.save(str(pdfPath))
        page.get_pixmap(dpi=pngDpi, alpha=False).save(str(pngPath))

    return (pdfPath, pngPath)

################################################################################
def displayPDFInNotebook(
    pdfPath : Path,
    dpi : int = 150, 
    width=900
    ) -> None :
    ''' 
    display a pdf in a jupyter notebook
    should work in vscode and jupyterlabs
    '''
    from IPython.display import Image as NotebookImage
    from IPython.display import display as notebookDisplay
    
    with fitz.open(pdfPath) as document:
        preview = document[0].get_pixmap(dpi=dpi)
        notebookDisplay(NotebookImage(data=preview.tobytes("png"), width=width))


try:
    from PIL import Image  # Pillow is commonly present with matplotlib
except ImportError:
    Image = None  # If missing, you can `pip install pillow`

PlotLike = Union[ConfusionMatrixDisplay, Figure, Axes]

################################################################################
def _fig_to_rgb_array_nondestructive(fig: Figure, dpi: Optional[int] = 150) -> np.ndarray:
    """
    Convert a Figure to an RGB array WITHOUT modifying the figure:
    - no canvas replacement
    - no DPI changes on the original
    Uses fig.savefig() to a BytesIO buffer.
    """
    if Image is None:
        raise RuntimeError("Pillow is required to rasterize Figures safely. `pip install pillow`")

    buf = io.BytesIO()
    # Keep bbox as None to avoid auto-tight-layout side effects on appearance
    logger.warning(f'AEDWIP 8/5/26 changed code to save png and pdf check file names and types')
    fig.savefig(buf, format="png", dpi=dpi, bbox_inches=None)
    fig.savefig(buf, format="pdf", dpi=dpi, bbox_inches=None)
    buf.seek(0)
    arr = np.asarray(Image.open(buf).convert("RGB"))
    return arr  # H x W x 3

################################################################################
def getSubplotWidthAndHeight(layout :dict[str, Panel],
                             panelId : str ) -> tuple[float, float]:
    ''' 
    todo
    
    returns (widthMM, heightMm)
    '''
    widthMM = layout[panelId].widthMm
    print(f'layout[{panelId}].widthMm : {widthMM}')
    heightMm = layout[panelId].heightMm
    print(f'layout[{panelId}].heightMm : {heightMm}')
    
    return(widthMM, heightMm)

################################################################################
def _choose_grid(n: int, ncols: Optional[int]) -> tuple[int, int]:
    if ncols is None:
        ncols = 3 if n >= 4 else (2 if n >= 2 else 1)
    nrows = (n + ncols - 1) // ncols
    return nrows, ncols

################################################################################
def make_panel(
    plots: List[PlotLike],
    ncols: Optional[int] = None,
    figsize: tuple[float, float] = (14, 9),
    mode: Literal["auto", "uniform", "mixed"] = "auto",
    raster_dpi: int = 300, #150,
    cm_shared_scale: bool = True,
    hide_colorbars: bool = True,
    titles: Optional[Sequence[str]] = None,
    displayAsPrecentage : bool = True
) -> Figure:
    """
    Build a grid from plot objects WITHOUT modifying originals.

    Modes:
      - "uniform": all same type (all ConfusionMatrixDisplay OR all Figure OR all Axes)
      - "mixed": different types/sizes
      - "auto": detect and choose

    Non-mutation guarantees:
      - Figures/Axes are rasterized via savefig→BytesIO (no canvas/DPI changes).
      - ConfusionMatrixDisplay is re-instantiated on target axes.
      
    displayAsPrecentage
        if True and ConfusionMatrixDisplay display value as a %
        
    figsize
        (width,height) default(14, 9)
    """
    if not plots:
        raise ValueError("No plots provided.")

    all_cmd = all(isinstance(p, ConfusionMatrixDisplay) for p in plots)
    all_fig = all(isinstance(p, Figure) for p in plots)
    all_ax  = all(isinstance(p, Axes) for p in plots)

    if mode == "auto":
        mode = "uniform" if (all_cmd or all_fig or all_ax) else "mixed"

    n = len(plots)
    nrows, ncols = _choose_grid(n, ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=figsize)
    axes = np.atleast_2d(axes)

    def _set_title(i: int, ax: Axes):
        if titles and i < len(titles):
            ax.set_title(titles[i])

    if mode == "uniform" and all_cmd:
        vmax = None
        if cm_shared_scale:
            vmax = max(int(np.asarray(p.confusion_matrix).max()) for p in plots)
        for i, disp in enumerate(plots):
            r, c = divmod(i, ncols)
            ax = axes[r, c]
            # Recreate a fresh ConfusionMatrixDisplay to avoid mutating the original object
            new_disp = ConfusionMatrixDisplay(
                confusion_matrix=np.asarray(disp.confusion_matrix),
                display_labels=getattr(disp, "display_labels", None),
            )
            im_kw = {}
            if cm_shared_scale and vmax is not None:
                im_kw = {"vmin": 0, "vmax": vmax}
                
            if displayAsPrecentage:
                new_disp.plot(ax=ax, colorbar=(not hide_colorbars), 
                            values_format=".0f", # <- forces fixed-point (100 instead of 1e+02)
                            im_kw=im_kw)
                for txt in ax.texts:
                    txt.set_text(txt.get_text() + "%")                
            else :
                new_disp.plot(ax=ax, colorbar=(not hide_colorbars), 
                    # values_format=".0f",         # <- forces fixed-point (100 instead of 1e+02)
                    im_kw=im_kw)
                
            _set_title(i, ax)
        # Hide extras
        for j in range(n, nrows * ncols):
            r, c = divmod(j, ncols)
            axes[r, c].axis("off")
        fig.tight_layout()
        return fig

    # Uniform Figures or Axes → rasterize each safely
    if mode == "uniform" and (all_fig or all_ax):
        for i, item in enumerate(plots):
            r, c = divmod(i, ncols)
            ax = axes[r, c]
            src_fig = item.figure if isinstance(item, Axes) else item
            arr = _fig_to_rgb_array_nondestructive(src_fig, dpi=raster_dpi)
            ax.imshow(arr)
            ax.axis("off")
            _set_title(i, ax)
        for j in range(n, nrows * ncols):
            r, c = divmod(j, ncols)
            axes[r, c].axis("off")
        fig.tight_layout()
        return fig

    # Mixed: re-create CMD; rasterize other types
    for i, item in enumerate(plots):
        r, c = divmod(i, ncols)
        ax = axes[r, c]
        if isinstance(item, ConfusionMatrixDisplay):
            new_disp = ConfusionMatrixDisplay(
                confusion_matrix=np.asarray(item.confusion_matrix),
                display_labels=getattr(item, "display_labels", None),
            )
            new_disp.plot(ax=ax, colorbar=(not hide_colorbars))
        elif isinstance(item, Axes):
            arr = _fig_to_rgb_array_nondestructive(item.figure, dpi=raster_dpi)
            ax.imshow(arr)
            ax.axis("off")
        elif isinstance(item, Figure):
            arr = _fig_to_rgb_array_nondestructive(item, dpi=raster_dpi)
            ax.imshow(arr)
            ax.axis("off")
        else:
            ax.text(0.5, 0.5, f"Unsupported: {type(item).__name__}",
                    ha="center", va="center", transform=ax.transAxes)
            ax.axis("off")
        _set_title(i, ax)

    for j in range(n, nrows * ncols):
        r, c = divmod(j, ncols)
        axes[r, c].axis("off")

    fig.tight_layout()
    return fig



################################################################################
def removeLegends(
        listOfFigs : list[plt.figure],
        panelWidthInMM : int,
        legendHeightMm : int,
        nCols : int,
    ) -> tuple[list[plt.figure], plt.figure]: 
    ''' 
    utility function for combineFigures()
    
    removes the legend from each figure. Create a figure we can use
    as an optional legend argument in combineFigures() using the first
    legend from the first figure in the list with a figure 
    
    returns 
        (deepcopy(listOfFigs), legendFig)
        
    '''
    # Capture the shared entries from the first axes with a legend.
    legendHandles = None
    legendLabels = None
    
    retListOfFigs = deepcopy(listOfFigs)


    for fig in retListOfFigs:
        for ax in fig.axes:
            legend = ax.get_legend()

            if legend is not None:
                if legendHandles is None:
                    legendHandles, _ = ax.get_legend_handles_labels()
                    # Preserve your displayed labels, including shortened names.
                    legendLabels = [
                        text.get_text() for text in legend.get_texts()
                    ]

                legend.remove()

    if legendHandles is None:
        raise ValueError("No axes legend was found in the supplied figures.")

    # Create a full-width legend strip at its final physical size.

    retLegendFig = plt.figure(
        figsize=(panelWidthInMM / 25.4, legendHeightMm / 25.4)
    )

    retLegendFig.legend(
        handles=legendHandles,
        labels=legendLabels,
        loc="center",
        ncol=nCols,
        fontsize=7,
        frameon=False,
    )   
    
    return (retListOfFigs, retLegendFig)


################################################################################  
def setFigureFont(fig, fontFamily="Arial"):
    ''' 
    make sure all figures use the correct final artwork compliant font
    '''
    for text in fig.findobj(match=Text):
        text.set_fontfamily(fontFamily)
        
################################################################################
# def exampleOfCombinedFiguresWithOptionalLegend():
    
#     panelWidthInMM = 180
#     panelHeightInMM = 100 # data grid only
#     legendHeightMm = 18
#     nPanelColumns = 3
#     nPanelRows = 2

#     originalFigList = createFigs(
#                                             data, 
#                                             panelWidthInMM=panelWidthInMM,
#                                             panelHeightInMM=panelHeightInMM,
#                                             nPanelColumns=nPanelColumns,
#                                             nPanelRows=nPanelRows,
#                                         )

    # figList, legendFig = removeLegends(originalFigList,
    #                       panelWidthInMM=F_r2_2_panelWidthInMM,
    #                       legendHeightMm=F_r2_2_legendHeightMm,
    #                       nCols=3, # 3 items in a single row F_r2_2_nPanelColumns
    #                       )
#     t = combineFigures(
#                         figList,
#                         F_r5_2OutputStem,
#                         widthMm=panelWidthInMM,
#                         heightMm=panelHeightInMM,
#                         nRows=nPanelRows,
#                         nCols=nPanelColumns,
#                         legendFigure=legendFig,
#                         legendHeightMm=legendHeightMm
#                         )

#     def plotBiotypeLegend(
#         biotypeColors,
#         panelWidthInMM=180,
#         panelHeightInMM=150,
#         nPanelColumns=2,
#         nPanelRows=3,
#         fontSize=7,
#         legendColumns=1,
#         ) -> tuple[plt.figure, plt.axes]:
#         """Create a legend-only subplot. Return (fig, ax)."""

#         figsize = (
#             (panelWidthInMM / nPanelColumns) / 25.4,
#             (panelHeightInMM / nPanelRows) / 25.4,
#         )

#         fig, ax = plt.subplots(figsize=figsize)
#         ax.set_axis_off()

#         # Preserve dictionary order and the exact supplied colors.
#         handles = [
#             Patch(
#                 facecolor=color,
#                 edgecolor="none",
#                 label=biotype.replace("_", " "),
#             )
#             for biotype, color in biotypeColors.items()
#         ]

#         ax.legend(
#             handles=handles,
#             loc="center",
#             bbox_to_anchor=(0.5, 0.5),
#             bbox_transform=fig.transFigure,
#             # ncol=1,
#             ncol=legendColumns,
#             frameon=False,
#             fontsize=fontSize,
#             handlelength=1.4,
#             handleheight=0.9,
#             handletextpad=0.7,
#             labelspacing=0.5,
#             borderaxespad=0,
#         )

#         return (fig, ax)



################################################################################
def validatePanelLayout(
    panels: dict[str, Panel],
    widthMm: float,
    heightMm: float,
    legendFigure=None,
    legendHeightMm: float = 10,
    maxWidthMm: float = 180,
    maxHeightMm: float = 185,
    requireFigures: bool = True,
) -> None:
    """
    Validate layout before assembly. Raise ValueError on failure.

    heightMm excludes the optional legend strip.
    maxHeightMm applies to the complete figure, including that strip.

    Checks layout and source canvas sizes, not full artwork compliance.
    Does not modify figures or panel positions.
    """
    errors = []
    validPanels = []
    tolerance = 1e-6

    if not panels:
        errors.append("No panels were supplied.")

    if not all(isfinite(v) and v > 0 for v in (widthMm, heightMm)):
        raise ValueError("Canvas width and height must be finite and positive.")

    totalHeightMm = heightMm

    if legendFigure is not None:
        if not isfinite(legendHeightMm) or legendHeightMm <= 0:
            errors.append("Legend height must be finite and positive.")
        else:
            totalHeightMm += legendHeightMm

    if widthMm > maxWidthMm + tolerance:
        errors.append(
            f"Canvas width {widthMm:g} mm exceeds {maxWidthMm:g} mm."
        )

    if totalHeightMm > maxHeightMm + tolerance:
        errors.append(
            f"Total height {totalHeightMm:g} mm exceeds "
            f"{maxHeightMm:g} mm, including the legend."
        )

    def checkFigureSize(figure, width, height, name):
        actualWidth, actualHeight = figure.get_size_inches() * 25.4

        if abs(actualWidth - width) > 0.1 or abs(actualHeight - height) > 0.1:
            errors.append(
                f"{name}: source canvas is "
                f"{actualWidth:.1f} × {actualHeight:.1f} mm; "
                f"placement is {width:g} × {height:g} mm. "
                "Regenerate the figure at its placement size to avoid scaling."
            )

    for key, panel in panels.items():
        name = f"Panel {key!r}"
        coordinates = (
            panel.leftMm,
            panel.topMm,
            panel.widthMm,
            panel.heightMm,
        )

        if not all(isfinite(v) for v in coordinates):
            errors.append(f"{name}: coordinates must be finite.")
            continue

        if panel.widthMm <= 0 or panel.heightMm <= 0:
            errors.append(f"{name}: width and height must be positive.")
            continue

        validPanels.append((key, panel))

        right = panel.leftMm + panel.widthMm
        bottom = panel.topMm + panel.heightMm

        if panel.leftMm < 0 or panel.topMm < 0:
            errors.append(f"{name}: left and top must be nonnegative.")

        if right > widthMm + tolerance:
            errors.append(
                f"{name}: right edge is {right:g} mm; "
                f"canvas width is {widthMm:g} mm."
            )

        if bottom > heightMm + tolerance:
            errors.append(
                f"{name}: bottom edge is {bottom:g} mm; "
                f"data-area height is {heightMm:g} mm."
            )

        if panel.figure is None:
            if requireFigures:
                errors.append(f"{name}: no figure assigned.")
        else:
            logger.debug(f'AEDWIP begin checkFigureSize() panel.label : {panel.label}')
            checkFigureSize(
                panel.figure, panel.widthMm, panel.heightMm, name
            )
            logger.debug(f'AEDWIP end checkFigureSize() panel.label : {panel.label}')


    # Touching edges are allowed; overlapping areas are not.
    for (keyA, a), (keyB, b) in combinations(validPanels, 2):
        overlapWidth = (
            min(a.leftMm + a.widthMm, b.leftMm + b.widthMm)
            - max(a.leftMm, b.leftMm)
        )
        overlapHeight = (
            min(a.topMm + a.heightMm, b.topMm + b.heightMm)
            - max(a.topMm, b.topMm)
        )

        if overlapWidth > tolerance and overlapHeight > tolerance:
            errors.append(
                f"Panels {keyA!r} and {keyB!r} overlap by "
                f"{overlapWidth:g} × {overlapHeight:g} mm."
            )

    if legendFigure is not None and legendHeightMm > 0:
        checkFigureSize(
            legendFigure, widthMm, legendHeightMm, "Legend"
        )

    requiredWidth = max(
        (p.leftMm + p.widthMm for _, p in validPanels), default=0
    )
    requiredHeight = max(
        (p.topMm + p.heightMm for _, p in validPanels), default=0
    )

    if errors:
        raise ValueError(
            "Invalid panel layout:\n- "
            + "\n- ".join(errors)
            + f"\nCurrent panel extents: "
              f"{requiredWidth:g} × {requiredHeight:g} mm."
            + "\nIncreasing the canvas alone does not resolve overlaps."
        )

    print(
        f"Layout valid: {widthMm:g} × {totalHeightMm:g} mm total; "
        f"panels occupy extents of "
        f"{requiredWidth:g} × {requiredHeight:g} mm."
    )