#!/usr/bin/env python3


def main():
    # Keep infer_bin.py as the stable entrypoint used by scripts/inference_wrapper.sh
    # but delegate the actual implementation to infer_llr.py.
    from infer_llr import main as _main
    return _main()


if __name__ == "__main__":
    raise SystemExit(main())
