# **ASTRAEA — COSMIC PULSE**

## **PRODUCT & FEATURE BLUEPRINT**

**NASA Space Apps Challenge 2026 · SPHEREx & Planet X · Advanced**

**Product tagline:** The Living History of the Universe

**Product positioning**

An AI-assisted interactive astronomical discovery platform that transforms SPHEREx's multi-epoch infrared observations into an accessible experience for detecting, comparing, and exploring cosmic objects.

# **1. PRODUCT OVERVIEW**

## **1.1. Core Idea**

ASTRAEA là nền tảng khám phá thiên văn tương tác, kết hợp dữ liệu quan sát đa thời điểm (Multi-Epoch Observations) từ sứ mệnh NASA SPHEREx với công nghệ xử lý dữ liệu, AI và trải nghiệm khám phá trực quan.

Nền tảng giúp người dùng:

- Khám phá các vùng trời thông qua tọa độ hoặc bản đồ thiên văn.

- Quan sát cùng một vùng trời qua nhiều lần quan sát SPHEREx.

- Phát hiện những thiên thể có dấu hiệu thay đổi vị trí hoặc độ sáng.

- So sánh hình ảnh giữa các thời điểm.

- Nhận diện các ứng viên dị thường cần được kiểm tra.

- Khám phá dữ liệu phổ hồng ngoại khi dữ liệu phù hợp.

- Kiểm tra bằng chứng khoa học và chất lượng dữ liệu.

- Lưu giữ và chia sẻ những khám phá.

ASTRAEA hướng đến cả người dùng phổ thông, người yêu thiên văn, sinh viên và nhà nghiên cứu.

## **1.2. Problem**

SPHEREx cung cấp lượng lớn dữ liệu quan sát hồng ngoại qua nhiều chu kỳ quét bầu trời. Việc nhận diện những thiên thể có chuyển động hoặc biến thiên trong khối lượng dữ liệu này đòi hỏi công cụ hiển thị, so sánh và phân tích phù hợp.

Các vấn đề cần giải quyết:

- Dữ liệu quan sát qua nhiều thời điểm khó so sánh thủ công.

- Những thay đổi nhỏ về vị trí hoặc độ sáng có thể khó nhận biết.

- Người dùng phổ thông khó tiếp cận dữ liệu thiên văn chuyên sâu.

- Việc xác định một tín hiệu bất thường có thực sự đáng quan tâm hay không cần thêm bằng chứng và thông tin chất lượng dữ liệu.

## **1.3. Solution**

ASTRAEA cung cấp một quy trình khám phá thống nhất:

**Discover → Retrieve → Observe → Animate → Compare → Detect → Investigate → Understand → Collect**

Người dùng chọn vùng trời hoặc đối tượng, hệ thống truy xuất dữ liệu SPHEREx, hiển thị chuỗi ảnh qua các lần quan sát, hỗ trợ so sánh và cung cấp thông tin để kiểm tra các ứng viên có dấu hiệu bất thường.

## **1.4. Core UX Flow**

1.  Người dùng chọn một đối tượng từ danh sách khám phá hoặc nhập RA/Dec.

2.  Hệ thống xác định vùng trời và truy xuất các quan sát SPHEREx có sẵn.

3.  Người dùng xem ảnh và danh sách các Epoch.

4.  Cosmic Motion Animation phát chuỗi ảnh theo thời gian.

5.  Người dùng dùng Time Slider, Blink hoặc Difference Mode để so sánh.

6.  Hệ thống hiển thị các nguồn sáng hoặc ứng viên có dấu hiệu thay đổi nếu pipeline xử lý hỗ trợ.

7.  Người dùng mở Technical Inspector để kiểm tra dữ liệu và bằng chứng.

8.  Người dùng lưu hoặc chia sẻ khám phá.

# **2. SYSTEM FEATURE ARCHITECTURE**

Hệ thống được chia thành ba nhóm:

- **Group I — Core MVP:** Các chức năng cốt lõi để đáp ứng mục tiêu NASA.

- **Group II — Advanced Intelligence:** Phân tích khoa học nâng cao và AI.

- **Group III — Immersive Experience:** Khả năng tiếp cận cộng đồng, kể chuyện và sáng tạo.

Priority:

- P0: Critical MVP.

- P1: Important extension.

- P2: Advanced / roadmap.

P0 thể hiện ưu tiên triển khai sản phẩm, không có nghĩa mọi mục trong Potential Considerations của NASA đều là yêu cầu bắt buộc.

# **GROUP I — CORE MVP**

## **1. Interactive Sky Map & Coordinate Explorer**

**Priority: P0**

Mục tiêu: Xác định vùng trời cần khám phá.

Chức năng:

- Interactive Sky Map.

- Pan & Zoom.

- Nhập Right Ascension (RA) và Declination (Dec).

- Tìm kiếm theo tên thiên thể.

- Chọn vùng trời trực tiếp trên bản đồ.

- Chuyển đổi tên thiên thể thành tọa độ khi có định danh phù hợp.

- Field of View (FOV) Controller.

- Hiển thị tọa độ trung tâm và phạm vi quan sát.

## **2. SPHEREx Multi-Epoch Image Viewer**

**Priority: P0**

Mục tiêu: Hiển thị các ảnh của cùng một vùng trời qua nhiều lần quan sát.

Chức năng:

- Truy xuất ảnh SPHEREx.

- Hiển thị các Epoch có sẵn.

- Observation Timeline.

- Time Slider.

- Smooth Cross-Fade.

- Auto-Play.

- Playback Speed.

- Hiển thị Observation Timestamp.

- Hỗ trợ nhiều Epoch khi dữ liệu có sẵn.

Yêu cầu xử lý:

- Kiểm tra hệ tọa độ.

- Kiểm tra kích thước và độ phân giải ảnh.

- Căn chỉnh ảnh khi cần so sánh.

- Chuẩn hóa cường độ khi phù hợp.

- Phân biệt dữ liệu gốc và ảnh đã xử lý.

## **3. Cosmic Motion Animation**

**Priority: P0 — Signature Core Feature**

Mục tiêu: Phát chuỗi ảnh của một vùng trời dưới dạng video ngắn hoặc animation liên tục, theo đúng mục tiêu của NASA.

Chức năng:

- Play / Pause.

- Loop.

- Playback Speed.

- Time Slider.

- Observation Timestamp.

- Continuous Playback.

- Blink Animation.

- Cross-Fade Animation.

Có thể bổ sung Export GIF/Video nếu đủ thời gian.

Nguyên tắc:

- Sử dụng các lần quan sát thực tế.

- Không trình bày các khung hình nội suy như quan sát thật.

- Phân biệt tốc độ phát với tốc độ chuyển động vật lý.

- Không tạo cảm giác có dữ liệu tại những Epoch không tồn tại.

## **4. Blink Comparator & Image Comparison**

**Priority: P0 for Blink; P1 for advanced processing**

Mục tiêu: Hỗ trợ người dùng nhận diện sự thay đổi giữa các lần quan sát.

| **Mode**        | **Function**                         |
|-----------------|--------------------------------------|
| Blink Mode      | Nhấp nháy giữa hai Epoch             |
| Side-by-side    | So sánh hai ảnh cạnh nhau            |
| Overlay Mode    | Chồng hai ảnh với opacity điều chỉnh |
| Cross-Fade      | Chuyển đổi mượt giữa ảnh             |
| Difference Mode | Hiển thị ảnh sai khác                |

Backend processing pipeline:

1.  Image Registration.

2.  Image Normalization.

3.  Image Subtraction.

4.  Source Detection.

5.  Candidate Identification.

6.  Candidate Tracking khi đủ dữ liệu.

Phân biệt hai nhóm:

- Moving Object Candidates.

- Variable Source Candidates.

Các khác biệt do nhiễu, sai số căn chỉnh, saturation và contamination phải được xem xét.

## **5. Metadata & Technical Inspector**

**Priority: P0**

Mục tiêu: Giúp người dùng hiểu nguồn gốc và chất lượng dữ liệu.

Thông tin hiển thị khi nguồn dữ liệu hỗ trợ:

- Object / Candidate ID.

- RA/Dec.

- Observation Timestamp.

- Epoch.

- Instrument / Detector.

- Data Product Type.

- Wavelength Band.

- Flux / Brightness.

- Signal-to-Noise Ratio.

- Contamination Flags.

- Quality Flags.

- Missing Data Indicators.

Phân loại:

- Known Object.

- Anomaly Candidate.

- Unclassified Source.

- Confirmed Object — chỉ khi được xác nhận bằng nguồn đáng tin cậy.

Có hai mức hiển thị:

- Public View: Nội dung dễ hiểu.

- Technical Inspector: Thông tin chuyên sâu.

## **6. NASA Benchmark Objects**

**Priority: P0**

Tạo khu vực NASA TEST OBJECTS.

| **Object**     | **Test Case**             |
|----------------|---------------------------|
| Barnard's Star | Proper Motion             |
| Epsilon Indi B | Brown Dwarf Binary System |
| WISE 0855      | Fast-moving Brown Dwarf   |
| RW Aurigae     | Variable Star             |
| R Draconis     | Variable Star             |

Mỗi đối tượng có:

- Tọa độ định danh.

- Nút mở vùng quan sát.

- Trạng thái dữ liệu.

- Các Epoch có sẵn.

- Kết quả so sánh nếu dữ liệu hỗ trợ.

- Metadata.

- Kết quả phân tích của hệ thống nếu triển khai được.

Không giả định cả năm đối tượng đều có cùng số Epoch hoặc đều có dữ liệu SPHEREx phù hợp. Backend cần xác minh từng trường hợp.

## **7. Guided Discovery & Curated Recommendation**

**Priority: P1**

Mục tiêu: Giúp người mới bắt đầu khám phá mà không cần biết tọa độ hoặc tên thiên thể.

Chức năng:

**Featured Discoveries**

- Danh sách đối tượng được tuyển chọn.

**Explore by Phenomenon**

- Moving Objects.

- Variable Stars.

- Brown Dwarfs.

- Asteroids and Comets.

- Unusual Sources.

**Discovery of the Day**

- Một đối tượng hoặc vùng trời được giới thiệu.

- Câu hỏi gợi mở.

- Nút Explore.

**Popular Searches**

- Chỉ sử dụng khi có dữ liệu lượt tìm kiếm thực tế.

- Nếu chưa có, dùng Featured Discoveries hoặc Recommended Targets.

**Personalized Recommendation — V2**

- Dựa trên lịch sử và sở thích khám phá.

V1 chỉ cần danh sách dữ liệu được tuyển chọn, không cần xây dựng machine learning recommendation.

## **8. User Guide**

**Priority: P1**

- How to Explore the Sky.

- How to Enter RA/Dec.

- How to Select FOV.

- How to Use Time Slider.

- How to Play Cosmic Motion Animation.

- How to Compare Images.

- How to Interpret Anomalies.

- How to Read Spectra.

- How to Understand Confidence and Data Quality.

# **GROUP II — ADVANCED INTELLIGENCE**

## **1. AI Anomaly Detection Engine**

**Priority: P1**

Mục tiêu: Đánh dấu những nguồn sáng hoặc đặc điểm dữ liệu khác biệt cần được kiểm tra.

Phương pháp:

- Statistical Outlier Detection.

- Unsupervised Anomaly Detection.

- Image Difference Analysis.

- Feature-based Source Detection.

Đầu ra:

- Candidate Coordinates.

- Anomaly Type.

- Anomaly Score.

- Signal-to-Noise Ratio.

- Observation Epoch.

- Difference Visualization.

- Data Quality Indicators.

Anomaly score không đồng nghĩa xác suất một thiên thể mới tồn tại.

## **2. Moving Object & Variable Source Analysis**

**Priority: P1**

**Moving Object Analysis**

- So sánh vị trí giữa các Epoch.

- Đánh giá positional change.

- Kiểm tra uncertainty.

- Theo dõi ứng viên khi đủ dữ liệu.

**Variable Source Analysis**

- Đo lường brightness variation.

- Tạo light curve khi đủ dữ liệu.

- Kiểm tra chất lượng phép đo.

- Phân biệt biến thiên thật với sai khác do thiết bị và xử lý.

## **3. Uncertainty-aware Candidate Ranking**

**Priority: P1/P2**

Tiêu chí:

- Anomaly Score.

- Positional Change.

- Brightness Variation.

- Signal-to-Noise Ratio.

- Data Quality.

- Measurement Uncertainty.

Đầu ra:

- Priority Candidate List.

- Review Status.

- Supporting Evidence.

- Uncertainty Indicators.

Không đồng nhất model confidence với xác suất một thiên thể mới thực sự tồn tại.

## **4. 102-Band Infrared Spectral Explorer**

**Priority: P1, conditional on data availability**

Chức năng:

- Hiển thị phổ theo bước sóng.

- Chọn các wavelength bands.

- So sánh cường độ giữa các dải.

- So sánh phổ giữa các Epoch khi phù hợp.

- Spectral Fingerprint.

- On-demand Spectrum Generation khi sản phẩm dữ liệu hỗ trợ.

Cần xác minh loại data product và metadata. Không mặc định mọi vùng trời đều có sẵn 102 ảnh độc lập tương ứng 102 dải.

## **5. Contrastive-driven Feature Topology**

**Priority: P2**

- Feature Extraction.

- Latent Representation.

- Similarity-based Clustering.

- Unusual Source Identification.

- UMAP Visualization.

Đây là hướng nghiên cứu nâng cao, không phải yêu cầu của MVP.

## **6. Sky Health Index**

**Priority: P2**

- Data Completeness.

- Image Quality.

- Contamination Level.

- Signal-to-Noise Quality.

- Valid Observation Coverage.

Các chỉ số cần công thức và dữ liệu hỗ trợ rõ ràng.

## **7. Multi-Survey Data Overlay**

**Priority: P1/P2**

Tích hợp Gaia, WISE và các catalog công khai khác.

Chức năng:

- Data Overlay.

- Cross-Matching.

- Position Comparison.

- Catalog Identification.

Phải xử lý khác biệt về tọa độ, epoch tham chiếu, độ phân giải và sai số vị trí.

## **8. AI Cosmic Guide**

**Priority: P2**

Trợ lý AI hỗ trợ giải thích:

- What am I looking at?

- What changed?

- Why does this object appear to move?

- What does this spectrum tell us?

- Why has this source been flagged?

AI phải trả lời dựa trên dữ liệu và kết quả phân tích thực tế, nêu rõ giới hạn khi thiếu bằng chứng.

# **GROUP III — IMMERSIVE EXPERIENCE**

## **1. Art-house Cosmic Storytelling**

**Priority: P1/P2**

- Cosmic Stories.

- Featured Celestial Objects.

- Guided Exploration.

- Historical Discovery Context.

- Thematic Discovery Journeys.

Nội dung phân biệt:

- Scientific Facts.

- Data-derived Results.

- Historical Context.

- Artistic Interpretation.

## **2. Cosmic Casebook & Discovery Journal**

**Priority: P1/P2**

**Cosmic Casebook:** Hồ sơ theo từng đối tượng hoặc vùng trời, gồm ảnh, tọa độ, Epoch, bằng chứng, metadata, kết quả và ghi chú.

**Discovery Journal:** Nhật ký theo hành trình khám phá, gồm vùng trời đã ghé thăm, câu hỏi, ứng viên đã kiểm tra và suy nghĩ cá nhân.

## **3. Cosmic Postcard & Static Image Export**

**Priority: P2**

Các dạng ảnh:

- Cosmic Portrait.

- Then & Now.

- Cosmic Difference.

- My Discovery.

Hỗ trợ xuất PNG nếu triển khai được.

Ảnh cần có thông tin nguồn, timestamp, tọa độ, chú giải màu và trạng thái xử lý.

Đây là đầu ra bổ sung, không thay thế Cosmic Motion Animation.

## **4. Cosmic Soundscape**

**Priority: P2 / Roadmap**

Chuyển đổi đặc trưng dữ liệu thành âm thanh thông qua Data Sonification.

Ví dụ:

- Signal Intensity → Amplitude.

- Wavelength → Pitch theo mapping xác định.

- Temporal Variation → Rhythm.

Không mô tả đây là âm thanh thực thu được từ vũ trụ.

## **5. Shareable Discovery**

**Priority: P1/P2**

- Chia sẻ vùng trời qua URL.

- Lưu tọa độ và FOV trong URL.

- Chia sẻ trạng thái quan sát nếu hỗ trợ.

- Xuất ảnh so sánh.

- Chia sẻ hồ sơ khám phá.

# **3. SIGNATURE EXPERIENCE — COSMIC PULSE**

Cosmic Pulse là trải nghiệm trung tâm của ASTRAEA.

**Discover → Observe → Travel Through Time → Compare → Detect → Investigate → Understand → Collect**

Đây là hành trình biến dữ liệu quan sát thành trải nghiệm khám phá.

Core interactions:

- Cosmic Motion Animation.

- Observation Timeline.

- Time Slider.

- Blink Comparator.

- Cross-Fade.

- Difference Mode.

- Technical Inspector.

# **4. SYSTEM ARCHITECTURE**

USER INTERFACE

- Interactive Sky Map.

- Coordinate Search.

- Featured Discoveries.

- Guided Discovery.

↓

COORDINATE & TARGET RESOLVER

- RA/Dec Validation.

- Object Name Resolution.

- FOV Selection.

↓

DATA RETRIEVAL LAYER

- SPHEREx Data Access.

- Observation Metadata.

- Multi-Epoch Images.

- Spectral Data when available.

↓

SCIENTIFIC PROCESSING LAYER

- FITS Processing.

- Image Registration.

- Image Normalization.

- Image Difference.

- Source Detection.

↓

SCIENTIFIC INTELLIGENCE LAYER

- Motion Analysis.

- Brightness Variation.

- Anomaly Detection.

- Candidate Ranking.

- Uncertainty Estimation.

↓

DISCOVERY EXPERIENCE LAYER

- Cosmic Motion Animation.

- Observation Viewer.

- Candidate Explorer.

- Technical Inspector.

- Cosmic Stories.

- Casebook.

## **Architectural Principles**

- Tách Data Retrieval, Processing, AI Analysis và Presentation.

- Backend chịu trách nhiệm xử lý dữ liệu khoa học.

- Frontend chịu trách nhiệm hiển thị và tương tác.

- Metadata và provenance cần được truyền cùng dữ liệu.

- Mock data phải được gắn nhãn rõ ràng.

- Phải có Loading, Empty, Error và Unavailable States.

## **Suggested Internal APIs**

Đây là endpoint nội bộ của ASTRAEA, không phải API chính thức của NASA:

- GET /api/search/coordinates

- GET /api/search/object

- GET /api/observations

- GET /api/observations/{id}

- GET /api/images/{id}

- POST /api/comparison

- GET /api/spectra/{id}

- GET /api/candidates

- GET /api/objects/{id}

Frontend và Backend cần thống nhất input/output JSON, đơn vị RA/Dec, timestamp, image format, Epoch ordering, quality flags và error responses.

# **5. NASA REQUIREMENT TRACEABILITY**

| **NASA objective / consideration** | **ASTRAEA feature**     | **Priority** |
|------------------------------------|-------------------------|--------------|
| RA/Dec Input                       | Coordinate Explorer     | P0           |
| SPHEREx Images                     | Multi-Epoch Viewer      | P0           |
| Video / Animation                  | Cosmic Motion Animation | P0           |
| Changes over Time                  | Timeline / Blink        | P0           |
| Moving Objects                     | Motion Analysis         | P1           |
| Brightness Changes                 | Variability Analysis    | P1           |
| Image Addition / Subtraction       | Image Arithmetic        | P1           |
| FOV                                | FOV Controller          | P0           |
| Five Benchmark Objects             | NASA Benchmark          | P0           |
| Metadata / Contamination Flags     | Technical Inspector     | P0           |
| New Data Integration               | Modular Pipeline        | P1           |
| Gaia / WISE                        | Multi-Survey Overlay    | P1/P2        |
| Spectra on Demand                  | Spectral Explorer       | P1/P2        |
| URL Coordinate Query               | Shareable Sky View      | P1           |
| AI Anomaly Detection               | Discovery Engine        | P1           |
| Contrastive Learning               | Feature Topology        | P2           |
| Storytelling / Audio               | Immersive Experience    | P2           |

# **6. 48-HOUR IMPLEMENTATION ROADMAP**

## **Hour 0–8: Data Validation & Foundation**

Backend + Astronomy:

- Kiểm tra SPHEREx access.

- Xác định data products.

- Kiểm tra Epoch, ảnh và metadata.

- Xác minh ít nhất một benchmark object.

- Chọn vùng trời demo có dữ liệu thật.

UI/UX + Frontend:

- Chốt Sitemap.

- Chốt User Flow.

- Thiết lập Design System.

- Thiết kế Sky Map và Viewer.

- Thống nhất API mock contract.

Deliverable: Một target có dữ liệu thực tế và giao diện nền tảng.

## **Hour 8–20: Core Astronomy Experience**

Backend:

- Coordinate Resolver.

- Observation Retrieval.

- Timestamp Extraction.

- FITS Processing.

- Image Validation.

- Image Registration nếu cần.

Frontend:

- Sky Map.

- Search.

- Observation Viewer.

- Timeline.

- Cosmic Motion Animation.

- Time Slider.

- Blink Comparator.

Deliverable: Một chuỗi quan sát thực tế có thể phát và so sánh.

## **Hour 20–34: Scientific Comparison & Discovery**

Backend:

- Image Normalization.

- Difference Image.

- Basic Source Detection nếu khả thi.

- Metadata API.

- Candidate Analysis nếu đủ dữ liệu.

Frontend:

- Difference Mode.

- Technical Inspector.

- Candidate Markers.

- Featured Discoveries.

- Guided Discovery.

- Loading / Error / Empty States.

Deliverable: Người dùng có thể quan sát, so sánh và kiểm tra bằng chứng.

## **Hour 34–48: Integration & Presentation**

- API integration.

- E2E testing.

- Benchmark testing.

- UI polish.

- Animation optimization.

- Demo preparation.

- Documentation and submission.

Optional:

- Basic Spectral Explorer.

- Static Image Export.

- AI Cosmic Guide.

- Advanced Anomaly Detection.

**Priority rule:** Không đánh đổi một quy trình khoa học thực sự hoạt động để lấy nhiều hiệu ứng hoặc module AI chưa được kiểm chứng.

# **7. MVP ACCEPTANCE CRITERIA**

| **Test**          | **Acceptance Criteria**                               |
|-------------------|-------------------------------------------------------|
| Coordinate Input  | Nhập RA/Dec và xác định đúng vùng trời                |
| FOV               | Thay đổi phạm vi quan sát                             |
| Data Retrieval    | Truy xuất dữ liệu SPHEREx hợp lệ                      |
| Multi-Epoch       | Hiển thị đúng Epoch có sẵn                            |
| Animation         | Phát chuỗi ảnh thành video/animation                  |
| Time Slider       | Chuyển được đến từng Epoch                            |
| Blink             | Luân phiên ảnh để nhận diện thay đổi                  |
| Comparison        | So sánh ảnh sau xử lý phù hợp                         |
| Metadata          | Hiển thị thông tin và chất lượng dữ liệu              |
| Benchmark         | Có quy trình kiểm thử năm đối tượng                   |
| Discovery Entry   | Người mới mở được đối tượng được tuyển chọn           |
| Guided Experience | Người mới thực hiện được một phiên khám phá           |
| Data Integrity    | Phân biệt dữ liệu thật, dữ liệu xử lý và dữ liệu demo |

## **Minimum Demo Success**

Một người dùng mới có thể:

1.  Mở ASTRAEA.

2.  Chọn Featured Discovery hoặc nhập RA/Dec.

3.  Mở vùng trời.

4.  Xem ảnh SPHEREx thực tế.

5.  Phát Cosmic Motion Animation.

6.  Dừng tại một Epoch.

7.  So sánh hai lần quan sát.

8.  Mở Technical Inspector.

9.  Hiểu điều gì đã thay đổi và giới hạn của dữ liệu.

# **8. PRODUCT ROADMAP**

## **V1 — The Living Sky Explorer**

- Coordinate Search.

- SPHEREx Retrieval.

- Multi-Epoch Viewer.

- Cosmic Motion Animation.

- Time Slider.

- Blink Comparator.

- Basic Image Comparison.

- Metadata Inspector.

- Featured Discoveries.

- Guided Discovery.

- Basic data-driven visual experience.

## **V2 — The Intelligent Cosmic Discovery Atlas**

- AI Anomaly Detection.

- Source Detection.

- Motion Tracking.

- Variability Analysis.

- Spectral Explorer.

- Cross-Matching.

- Candidate Ranking.

- Uncertainty Analysis.

- Contrastive Feature Topology.

- AI Cosmic Guide.

- Sky Health.

- Cosmic Casebook.

- Discovery Journal.

- Personalized Recommendations.

## **V3 — The Immersive Living Universe**

- Generative Cosmic Art.

- Cosmic Soundscape.

- Personal Cosmic Narratives.

- Advanced Storytelling.

- Immersive Projection.

- Community Discovery Network.

# **9. FINAL PRODUCT POSITIONING**

## **ASTRAEA — COSMIC PULSE**

### **The Living History of the Universe**

**Turning the invisible motion of the universe into an experience everyone can explore.**

Three Core Values:

**Scientific Discovery:** Khám phá những thay đổi của thiên thể qua dữ liệu quan sát thực tế.

**Artificial Intelligence:** Hỗ trợ nhận diện, phân tích và ưu tiên những ứng viên đáng được nghiên cứu thêm.

**Human-centered Exploration:** Biến dữ liệu thiên văn phức tạp thành một hành trình khám phá trực quan và dễ tiếp cận.

**Core UX Flow**

Discover → Retrieve → Animate → Compare → Detect → Investigate → Understand → Collect.
