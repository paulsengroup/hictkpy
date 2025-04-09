# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT


def _get_package_version(name: str, missing_ok: bool):
    try:
        from importlib.metadata import version

        return version(name)
    except:  # noqa
        if missing_ok:
            return None
        raise


from ._hictkpy import (
    Bin,
    BinTable,
    File,
    MultiResFile,
    PixelSelector,
    __doc__,
    __hictk_version__,
    cooler,
    hic,
    is_cooler,
    is_hic,
    is_mcool_file,
    is_scool_file,
)

__version__ = _get_package_version("hictkpy", missing_ok=False)
__all__ = [
    "__doc__",
    "Bin",
    "BinTable",
    "File",
    "MultiResFile",
    "PixelSelector",
    "is_cooler",
    "is_mcool_file",
    "is_scool_file",
    "is_hic",
    "cooler",
    "hic",
    "__hictk_version__",
]


def _try_report_successful_import():
    try:
        import os

        if "HICTK_NO_TELEMETRY" in os.environ or "HICTKPY_NO_TELEMETRY" in os.environ:
            return

        import platform

        from opentelemetry import trace
        from opentelemetry.exporter.otlp.proto.http.trace_exporter import (
            OTLPSpanExporter,
        )
        from opentelemetry.sdk.resources import (
            HOST_ARCH,
            OS_TYPE,
            SERVICE_NAME,
            SERVICE_VERSION,
            Resource,
        )
        from opentelemetry.sdk.trace import TracerProvider
        from opentelemetry.sdk.trace.export import (
            BatchSpanProcessor,
            ConsoleSpanExporter,
        )

        resource = Resource(
            attributes={
                SERVICE_NAME: "hictkpy",
                SERVICE_VERSION: __version__,
                OS_TYPE: platform.system(),
                HOST_ARCH: platform.machine(),
                "python.version": platform.python_version(),
                "pandas.version": _get_package_version("pandas", missing_ok=True),
                "numpy.version": _get_package_version("numpy", missing_ok=True),
                "pyarrow.version": _get_package_version("pyarrow", missing_ok=True),
                "scipy.version": _get_package_version("scipy", missing_ok=True),
                "schema": 1,
            }
        )

        provider = TracerProvider(resource=resource)
        processor = BatchSpanProcessor(
            OTLPSpanExporter(
                # endpoint="https://hictkpy-telemetry.paulsenlab.com:4318",
                endpoint="localhost:9",
            )
            # ConsoleSpanExporter()
        )
        provider.add_span_processor(processor)

        # Sets the global default tracer provider
        trace.set_tracer_provider(provider)

        # Creates a tracer from the global tracer provider
        tracer = trace.get_tracer("hictkpy")

        with tracer.start_as_current_span("module-imported") as span:
            span.set_status(trace.status.StatusCode.OK)

    except:  # noqa
        pass


_try_report_successful_import()
