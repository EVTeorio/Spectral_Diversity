Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

$root = Split-Path -Parent $PSScriptRoot
$outDir = Join-Path $root "Documents/Tables and Figures"
$outPath = Join-Path $outDir "19_alpha_hull_vs_convex_hull_concept.png"
New-Item -ItemType Directory -Force -Path $outDir | Out-Null

Add-Type -AssemblyName System.Drawing

function Font([single]$size, $style = [System.Drawing.FontStyle]::Regular) {
  [System.Drawing.Font]::new("Arial", $size, $style)
}

function Brush($color) {
  [System.Drawing.SolidBrush]::new($color)
}

function PenC($color, [single]$width = 2) {
  [System.Drawing.Pen]::new($color, $width)
}

function Draw-Text($g, [string]$text, [single]$x, [single]$y, [single]$size, $color, $style = [System.Drawing.FontStyle]::Regular) {
  $font = Font $size $style
  $brush = Brush $color
  $g.DrawString($text, $font, $brush, $x, $y)
  $font.Dispose()
  $brush.Dispose()
}

function Draw-Centered($g, [string]$text, [single]$x, [single]$y, [single]$w, [single]$h, [single]$size, $color, $style = [System.Drawing.FontStyle]::Regular) {
  $font = Font $size $style
  $brush = Brush $color
  $fmt = [System.Drawing.StringFormat]::new()
  $fmt.Alignment = [System.Drawing.StringAlignment]::Center
  $fmt.LineAlignment = [System.Drawing.StringAlignment]::Center
  $g.DrawString($text, $font, $brush, [System.Drawing.RectangleF]::new($x, $y, $w, $h), $fmt)
  $fmt.Dispose()
  $font.Dispose()
  $brush.Dispose()
}

function New-PointF([double]$x, [double]$y) {
  [System.Drawing.PointF]::new([single]$x, [single]$y)
}

function Shift-Points($points, [double]$dx, [double]$dy) {
  $points | ForEach-Object { New-PointF ($_.X + $dx) ($_.Y + $dy) }
}

function Draw-Panel($g, [string]$title, [string]$subtitle, [double]$x0, [bool]$alphaHull) {
  $ink = [System.Drawing.Color]::FromArgb(36, 44, 52)
  $muted = [System.Drawing.Color]::FromArgb(92, 100, 110)
  $axis = [System.Drawing.Color]::FromArgb(130, 140, 150)
  $pointFill = [System.Drawing.Color]::FromArgb(42, 86, 145)
  $convexFill = [System.Drawing.Color]::FromArgb(62, 70, 126, 185)
  $convexEdge = [System.Drawing.Color]::FromArgb(70, 126, 185)
  $alphaFill = [System.Drawing.Color]::FromArgb(70, 36, 150, 100)
  $alphaEdge = [System.Drawing.Color]::FromArgb(36, 150, 100)
  $emptyFill = [System.Drawing.Color]::FromArgb(85, 210, 170, 70)

  Draw-Centered $g $title $x0 165 900 45 32 $ink ([System.Drawing.FontStyle]::Bold)
  Draw-Centered $g $subtitle $x0 212 900 70 23 $muted

  $plotX = $x0 + 80
  $plotY = 315
  $plotW = 740
  $plotH = 600
  $g.DrawLine((PenC $axis 3), $plotX, $plotY + $plotH, $plotX + $plotW, $plotY + $plotH)
  $g.DrawLine((PenC $axis 3), $plotX, $plotY, $plotX, $plotY + $plotH)
  Draw-Centered $g "PC1" ($plotX + 300) ($plotY + $plotH + 38) 140 36 22 $muted

  $state = $g.Save()
  $g.TranslateTransform($plotX - 58, $plotY + 380)
  $g.RotateTransform(-90)
  Draw-Centered $g "PC2" 0 0 140 36 22 $muted
  $g.Restore($state)

  $points = @(
    (New-PointF 180 525), (New-PointF 225 480), (New-PointF 270 438), (New-PointF 325 408),
    (New-PointF 390 396), (New-PointF 455 408), (New-PointF 510 438), (New-PointF 560 482),
    (New-PointF 600 532), (New-PointF 532 552), (New-PointF 465 562), (New-PointF 398 558),
    (New-PointF 335 540), (New-PointF 280 510), (New-PointF 255 340), (New-PointF 320 300),
    (New-PointF 395 284), (New-PointF 470 300), (New-PointF 535 342), (New-PointF 505 372),
    (New-PointF 440 350), (New-PointF 372 350), (New-PointF 312 372)
  )
  $drawPoints = @(Shift-Points $points $plotX $plotY)

  $convex = @(
    (New-PointF 180 525), (New-PointF 255 340), (New-PointF 320 300), (New-PointF 395 284),
    (New-PointF 470 300), (New-PointF 535 342), (New-PointF 600 532), (New-PointF 532 552),
    (New-PointF 465 562), (New-PointF 398 558), (New-PointF 335 540), (New-PointF 280 510)
  )
  $alpha = @(
    (New-PointF 180 525), (New-PointF 225 480), (New-PointF 270 438), (New-PointF 325 408),
    (New-PointF 255 340), (New-PointF 320 300), (New-PointF 395 284), (New-PointF 470 300),
    (New-PointF 535 342), (New-PointF 510 438), (New-PointF 560 482), (New-PointF 600 532),
    (New-PointF 532 552), (New-PointF 465 562), (New-PointF 398 558), (New-PointF 335 540),
    (New-PointF 280 510)
  )

  if ($alphaHull) {
    $poly = @(Shift-Points $alpha $plotX $plotY)
    $g.FillPolygon((Brush $alphaFill), [System.Drawing.PointF[]]$poly)
    $g.DrawPolygon((PenC $alphaEdge 6), [System.Drawing.PointF[]]$poly)
    Draw-Centered $g "less empty spectral space" ($plotX + 180) ($plotY + 170) 320 44 21 $alphaEdge ([System.Drawing.FontStyle]::Bold)
  } else {
    $poly = @(Shift-Points $convex $plotX $plotY)
    $g.FillPolygon((Brush $convexFill), [System.Drawing.PointF[]]$poly)
    $g.DrawPolygon((PenC $convexEdge 6), [System.Drawing.PointF[]]$poly)
    $empty = @(
      (New-PointF 325 408), (New-PointF 390 396), (New-PointF 455 408),
      (New-PointF 510 438), (New-PointF 535 342), (New-PointF 470 300),
      (New-PointF 395 284), (New-PointF 320 300), (New-PointF 255 340)
    )
    $emptyShift = @(Shift-Points $empty $plotX $plotY)
    $g.FillPolygon((Brush $emptyFill), [System.Drawing.PointF[]]$emptyShift)
    Draw-Centered $g "empty spectral space`nincluded" ($plotX + 312) ($plotY + 318) 260 72 22 ([System.Drawing.Color]::FromArgb(150, 100, 25)) ([System.Drawing.FontStyle]::Bold)
  }

  foreach ($p in $drawPoints) {
    $g.FillEllipse((Brush ([System.Drawing.Color]::White)), $p.X - 8, $p.Y - 8, 16, 16)
    $g.FillEllipse((Brush $pointFill), $p.X - 6, $p.Y - 6, 12, 12)
  }
}

$w = 2200
$h = 1200
$bmp = [System.Drawing.Bitmap]::new($w, $h, [System.Drawing.Imaging.PixelFormat]::Format32bppArgb)
$g = [System.Drawing.Graphics]::FromImage($bmp)
$g.Clear([System.Drawing.Color]::Transparent)
$g.SmoothingMode = [System.Drawing.Drawing2D.SmoothingMode]::AntiAlias
$g.TextRenderingHint = [System.Drawing.Text.TextRenderingHint]::ClearTypeGridFit

$ink = [System.Drawing.Color]::FromArgb(36, 44, 52)
$muted = [System.Drawing.Color]::FromArgb(92, 100, 110)

Draw-Text $g "Alpha-hull area versus convex-hull area in PCA spectral space" 80 42 38 $ink ([System.Drawing.FontStyle]::Bold)
Draw-Text $g "Both use the same pixel scores; they differ in how tightly the boundary follows occupied spectral space." 80 96 24 $muted

Draw-Panel $g "Convex hull" "Smallest convex polygon around all pixels" 95 $false
Draw-Panel $g "Alpha hull" "Boundary can follow concave structure in the point cloud" 1200 $true

Draw-Text $g "In this study, alpha-hull area was used to describe the occupied extent of pixel spectra in vector-normalized PCA space." 120 1072 24 $muted

$bmp.Save($outPath, [System.Drawing.Imaging.ImageFormat]::Png)
$g.Dispose()
$bmp.Dispose()

Write-Host "Created $outPath"
