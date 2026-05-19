from __future__ import annotations

from typing import Optional, Sequence

from .config import ANCOMConfig, Subset
from .data import ANCOMData
from .heatmap import HeatmapPlotter
from .trajectory import TrajectoryPlotter


class ANCOMWorkflow:
    def __init__(self, config: ANCOMConfig):
        self.config = config
        self.data = ANCOMData(config)
        self.heatmap_plotter = HeatmapPlotter(self.data)
        self.trajectory_plotter = TrajectoryPlotter(self.data)

    def check(self, subset: Subset):
        """
        Show which table and ANCOM paths would be used.
        Useful before running plots on a new dataset.
        """
        return self.data.explain_paths(subset)

    def heatmap(self, subset: Subset, **kwargs):
        return self.heatmap_plotter.plot(subset=subset, **kwargs)

    def heatmaps(self, subsets: Sequence[Subset], **kwargs):
        results = []
        for subset in subsets:
            results.append(self.heatmap(subset, **kwargs))
        return results

    def trajectory(
        self,
        taxon_query: str,
        subset: Subset,
        comparison_levels: Optional[Sequence[str]] = None,
    ):
        return self.trajectory_plotter.line(
            taxon_query=taxon_query,
            subset=subset,
            comparison_levels=comparison_levels,
        )

    def trajectories(
        self,
        subset: Subset,
        taxa: Optional[Sequence[str]] = None,
        comparison_levels: Optional[Sequence[str]] = None,
    ):
        taxa = list(taxa) if taxa is not None else list(self.config.taxa)

        if not taxa:
            raise ValueError(
                "No taxa configured. Add taxa=[...] to ANCOMConfig "
                "or pass taxa=[...] to workflow.trajectories()."
            )

        results = []
        for taxon in taxa:
            results.append(
                self.trajectory(
                    taxon_query=taxon,
                    subset=subset,
                    comparison_levels=comparison_levels,
                )
            )
        return results

    def boxplot(
        self,
        taxon_query: str,
        subset: Subset,
        comparison_levels: Optional[Sequence[str]] = None,
        show_trend: bool = True,
        trend_order: int = 2,
    ):
        return self.trajectory_plotter.boxplot(
            taxon_query=taxon_query,
            subset=subset,
            comparison_levels=comparison_levels,
            show_trend=show_trend,
            trend_order=trend_order,
        )

    def boxplots(
        self,
        subset: Subset,
        taxa: Optional[Sequence[str]] = None,
        comparison_levels: Optional[Sequence[str]] = None,
        show_trend: bool = True,
        trend_order: int = 2,
    ):
        taxa = list(taxa) if taxa is not None else list(self.config.taxa)

        if not taxa:
            raise ValueError(
                "No taxa configured. Add taxa=[...] to ANCOMConfig "
                "or pass taxa=[...] to workflow.boxplots()."
            )

        results = []
        for taxon in taxa:
            results.append(
                self.boxplot(
                    taxon_query=taxon,
                    subset=subset,
                    comparison_levels=comparison_levels,
                    show_trend=show_trend,
                    trend_order=trend_order,
                )
            )
        return results
