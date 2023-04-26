from pathlib import Path
from typing import Optional

from pydantic import BaseModel, Field


class PostProcessing(BaseModel):
    sort_order: str = Field(...)


class ExomiserDB(BaseModel):
    exomiser_database: Optional[Path] = Field(None)
    exomiser_hg19_database: Optional[Path] = Field(None)
    exomiser_hg38_database: Optional[Path] = Field(None)


class LIRICALToolSpecificConfigurations(BaseModel):
    mode: str = Field(...)
    lirical_jar_executable: Path = Field(...)
    exomiser_db_configurations: ExomiserDB = Field(...)
    post_process: PostProcessing = Field(...)
