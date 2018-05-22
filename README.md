# DB設計

## usersテーブル

|Column|Type|Options|
|------|----|-------|
|name|string|null: false, index: true|
|email|string|null: false, unique: true|

### Association
- has_many :posts, dependent: :destroy
- has_many :likes, dependent: :destroy
- has_many :comments, dependent: :destroy
- has_many :tags, through: :user_tags
- has_many :user_tags, dependent: :destroy
- has_many :followed_users, through: :relationships, source: :followed_id
- has_many :relationships, foreign_key: "follower_id", dependent: :destroy
- has_many :reverse_relationships, foreign_key: "followed_id", clas_name: "Relationship", dependent: :destroy
- has_many :followers, through: :reverse_relationships, source: :follower_id
- has_many :groups, through: groups_users, dependent: :destroy
- has_many :group_users, dependent: :destroy


## postsテーブル

|Column|Type|Options|
|------|----|-------|
|content|text|null: false|
|image_name|string||
|user_id|integer|null :false, index: true, foreign_key: true|

### Association
- belongs_to :user
- has_one    :group
- has_many   :likes, dependent: :destroy
- has_many   :comments
- has_many   :tags, through: :post_tags
- has_many   :post_tags, dependent: :destroy

## groupsテーブル

|Column|Type|Options|
|------|----|-------|
|user_id|integer|null :false, index: true, foreign_key: true|
|post_id|integer|null :false, index: true, foreign_key: true|

### Association
- has_many :users, through: group_users
- has_many :group_messages
- has_many :group_users


## group_usersテーブル

|Column|Type|Options|
|------|----|-------|
|user_id|integer|null :false, index: true, foreign_key: true|
|post_id|integer|null :false, index: true, foreign_key: true|

### Association
- belongs_to :user
- belongs_to :group


## group_messagesテーブル

|Column|Type|Options|
|------|----|-------|
|content|text|null :false|
|image_name|string||
|group_id|integer|null :false, index: true, foreign_key: true|
|user_id|integer|null :false, index: true, foreign_key: true|

### Association
- belongs_to :group
- belongs_to :user


## likesテーブル

|Column|Type|Options|
|------|----|-------|
|user_id|integer|null :false, index: true, foreign_key: true|
|post_id|integer|null :false, index: true, foreign_key: true|

### Association
- belongs_to :user
- belongs_to :post


## commentsテーブル

|Column|Type|Options|
|------|----|-------|
|content|text|null :false|
|user_id|integer|null :false, index: true, foreign_key: true|
|post_id|integer|null :false, index: true, foreign_key: true|

### Association
- belongs_to :user
- belongs_to :post


## post_tagsテーブル

|Column|Type|Options|
|------|----|-------|
|tag_id|integer|null :false, index: true, foreign_key: true|
|post_id|integer|null :false, index: true, foreign_key: true|

### Association
- belongs_to :post
- belongs_to :tags

## user_tagsテーブル

|Column|Type|Options|
|------|----|-------|
|tag_id|integer|null :false, index: true, foreign_key: true|
|user_id|integer|null :false, index: true, foreign_key: true|

### Association
- belongs_to :user
- belongs_to :tag


## tagsテーブル

|Column|Type|Options|
|------|----|-------|
|name|string|null :false, index: true, unique: true|

### Association
- has_many :posts, through: :post_tags
- has_many :post_tags, dependent: :destroy
- has_many :users, through: :user_tags
- has_many :user_tags, dependent: :destroy


## relationshipsテーブル

|Column|Type|Options|
|------|----|-------|
|follower_id|integer|null :false, index: true, foreign_key: true, unique: true|
|followed_id|integer|null :false, index: true, foreign_key: true, unique: true|

### Association
- belongs_to :followed, class_name: "User"
- belongs_to :follower, class_name: "User"
















