-- Convenience views over the parquet store for interactive duckdb sessions:
--   duckdb -init /work/submodules/g4cal/scripts/views.sql
CREATE OR REPLACE VIEW runs AS
  SELECT * FROM read_parquet('/data/ilreco-work/store/runs/*.parquet');
CREATE OR REPLACE VIEW cal_hits AS
  SELECT * FROM read_parquet('/data/ilreco-work/store/cal_hits/*.parquet');
CREATE OR REPLACE VIEW cal_events AS
  SELECT * FROM read_parquet('/data/ilreco-work/store/cal_events/*.parquet');
-- events joined with their run metadata, calibration runs excluded
CREATE OR REPLACE VIEW ev AS
  SELECT e.*, r.preset, r.particle, r.e_min, r.e_max, r.gun_mode
  FROM cal_events e JOIN runs r USING (run_id)
  WHERE r.run_id NOT LIKE 'calib-%' AND r.run_id NOT LIKE 'disp-%';
