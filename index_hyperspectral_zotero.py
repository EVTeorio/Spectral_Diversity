from __future__ import annotations

import csv
import json
import os
import sqlite3
from datetime import datetime
from pathlib import Path


COLLECTION_NAME = "Hyperspectral Paper"


def one(cur: sqlite3.Cursor, sql: str, params: tuple = ()):
    rows = cur.execute(sql, params).fetchall()
    if len(rows) != 1:
        return None, rows
    return rows[0], rows


def table_columns(cur: sqlite3.Cursor, table: str) -> list[str]:
    return [row[1] for row in cur.execute(f"pragma table_info({table})").fetchall()]


def item_data(cur: sqlite3.Cursor, item_id: int) -> dict[str, str]:
    rows = cur.execute(
        """
        select f.fieldName, v.value
        from itemData d
        join fields f on f.fieldID = d.fieldID
        join itemDataValues v on v.valueID = d.valueID
        where d.itemID = ?
        """,
        (item_id,),
    ).fetchall()
    return {field: value for field, value in rows}


def creators(cur: sqlite3.Cursor, item_id: int) -> str:
    rows = cur.execute(
        """
        select c.lastName, c.firstName
        from itemCreators ic
        join creators c on c.creatorID = ic.creatorID
        where ic.itemID = ?
        order by ic.orderIndex
        """,
        (item_id,),
    ).fetchall()
    names = []
    for last, first in rows:
        if first and last:
            names.append(f"{last}, {first}")
        else:
            names.append(last or first or "")
    return "; ".join(name for name in names if name)


def main() -> None:
    zotero_dir = Path(os.environ["USERPROFILE"]) / "Zotero"
    db_path = zotero_dir / "zotero.sqlite"
    storage_dir = zotero_dir / "storage"
    out_dir = Path.cwd() / "zotero_index"
    out_dir.mkdir(exist_ok=True)

    con = sqlite3.connect(f"file:{db_path}?mode=ro&immutable=1", uri=True)
    cur = con.cursor()

    match, matches = one(
        cur,
        "select collectionID, collectionName, parentCollectionID, key from collections where collectionName = ?",
        (COLLECTION_NAME,),
    )
    if match is None:
        print(f"Could not find exactly one collection named {COLLECTION_NAME!r}. Matches:")
        for row in matches:
            print(row)
        like = cur.execute(
            "select collectionID, collectionName, parentCollectionID, key from collections where lower(collectionName) like ? order by collectionName",
            (f"%{COLLECTION_NAME.lower()}%",),
        ).fetchall()
        for row in like:
            print(row)
        raise SystemExit(1)

    collection_id, collection_name, parent_collection_id, collection_key = match

    direct_item_ids = [
        row[0]
        for row in cur.execute(
            "select itemID from collectionItems where collectionID = ?",
            (collection_id,),
        ).fetchall()
    ]

    item_ids = set(direct_item_ids)
    child_item_ids = [
        row[0]
        for row in cur.execute(
            f"select itemID from itemAttachments where parentItemID in ({','.join('?' for _ in direct_item_ids)})",
            tuple(direct_item_ids),
        ).fetchall()
    ] if direct_item_ids else []
    item_ids.update(child_item_ids)

    item_rows = []
    file_rows = []
    for item_id in sorted(item_ids):
        item_key, item_type = cur.execute(
            """
            select i.key, t.typeName
            from items i
            join itemTypes t on t.itemTypeID = i.itemTypeID
            where i.itemID = ?
            """,
            (item_id,),
        ).fetchone()
        data = item_data(cur, item_id)
        attachment = None
        if item_id in child_item_ids or item_id in direct_item_ids:
            attachment = cur.execute(
                "select parentItemID, linkMode, contentType, path from itemAttachments where itemID = ?",
                (item_id,),
            ).fetchone()

        parent_id = attachment[0] if attachment else None
        parent_data = item_data(cur, parent_id) if parent_id else {}
        parent_key = None
        parent_type = None
        parent_creators = ""
        if parent_id:
            parent_key, parent_type = cur.execute(
                """
                select i.key, t.typeName
                from items i
                join itemTypes t on t.itemTypeID = i.itemTypeID
                where i.itemID = ?
                """,
                (parent_id,),
            ).fetchone()
            parent_creators = creators(cur, parent_id)

        file_path_obj = None
        file_path = ""
        exists = ""
        size = ""
        modified = ""
        if attachment:
            raw_path = attachment[3] or ""
            if raw_path.startswith("storage:"):
                file_path_obj = storage_dir / item_key / raw_path.replace("storage:", "", 1)
            else:
                file_path_obj = Path(raw_path)
            file_path = str(file_path_obj)
            exists = str(file_path_obj.exists())
            if file_path_obj.exists():
                stat = file_path_obj.stat()
                size = str(stat.st_size)
                modified = datetime.fromtimestamp(stat.st_mtime).isoformat(timespec="seconds")

        base = {
            "collection": collection_name,
            "collection_id": collection_id,
            "collection_key": collection_key,
            "item_id": item_id,
            "item_key": item_key,
            "item_type": item_type,
            "title": parent_data.get("title") or data.get("title") or "",
            "attachment_title": data.get("title") if attachment else "",
            "date": data.get("date") or parent_data.get("date") or "",
            "creators": creators(cur, item_id) or parent_creators,
            "parent_item_id": parent_id or "",
            "parent_key": parent_key or "",
            "parent_type": parent_type or "",
            "content_type": attachment[2] if attachment else "",
            "attachment_path": attachment[3] if attachment else "",
            "resolved_file_path": file_path,
            "file_exists": exists,
            "file_size_bytes": size,
            "file_modified": modified,
        }
        item_rows.append(base)

        if attachment and file_path_obj is not None:
            folder = file_path_obj.parent if file_path_obj.name else file_path_obj
            if folder.exists() and folder.is_dir():
                for actual_file in sorted(folder.iterdir(), key=lambda p: p.name.lower()):
                    if not actual_file.is_file():
                        continue
                    stat = actual_file.stat()
                    file_rows.append(
                        {
                            "collection": collection_name,
                            "title": base["title"],
                            "attachment_title": base["attachment_title"],
                            "date": base["date"],
                            "creators": base["creators"],
                            "attachment_item_id": item_id,
                            "attachment_item_key": item_key,
                            "attachment_content_type": attachment[2] or "",
                            "storage_folder": str(folder),
                            "file_name": actual_file.name,
                            "extension": actual_file.suffix,
                            "file_path": str(actual_file),
                            "file_size_bytes": stat.st_size,
                            "file_modified": datetime.fromtimestamp(stat.st_mtime).isoformat(timespec="seconds"),
                            "is_primary_attachment": str(actual_file == file_path_obj),
                        }
                    )

    items_csv_path = out_dir / "hyperspectral_paper_zotero_items_and_attachments.csv"
    items_json_path = out_dir / "hyperspectral_paper_zotero_items_and_attachments.json"
    with items_csv_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(item_rows[0].keys()) if item_rows else [])
        writer.writeheader()
        writer.writerows(item_rows)
    items_json_path.write_text(json.dumps(item_rows, indent=2, ensure_ascii=False), encoding="utf-8")

    files_csv_path = out_dir / "hyperspectral_paper_zotero_files_only.csv"
    files_json_path = out_dir / "hyperspectral_paper_zotero_files_only.json"
    with files_csv_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(file_rows[0].keys()) if file_rows else [])
        writer.writeheader()
        writer.writerows(file_rows)
    files_json_path.write_text(json.dumps(file_rows, indent=2, ensure_ascii=False), encoding="utf-8")

    print(f"collection_id={collection_id}")
    print(f"direct_items={len(direct_item_ids)}")
    print(f"child_attachments={len(child_item_ids)}")
    print(f"item_and_attachment_rows={len(item_rows)}")
    print(f"files_only_rows={len(file_rows)}")
    print(f"files_csv={files_csv_path}")
    print(f"files_json={files_json_path}")
    print(f"items_csv={items_csv_path}")
    print(f"items_json={items_json_path}")


if __name__ == "__main__":
    main()
